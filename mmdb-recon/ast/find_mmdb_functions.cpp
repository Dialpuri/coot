#include "clang/AST/AST.h"
#include "clang/AST/ASTConsumer.h"
#include "clang/AST/RecursiveASTVisitor.h"
#include "clang/Frontend/ASTConsumers.h"
#include "clang/Frontend/CompilerInstance.h"
#include "clang/Frontend/FrontendActions.h"
#include "clang/Lex/PPCallbacks.h"
#include "clang/Lex/Preprocessor.h"
#include "clang/Tooling/ArgumentsAdjusters.h"
#include "clang/Tooling/CommonOptionsParser.h"
#include "clang/Tooling/Tooling.h"
#include "llvm/Support/CommandLine.h"
#include "llvm/Support/JSON.h"

#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <stdarg.h>

using namespace clang;
using namespace clang::tooling;

// Structure to hold information about mmdb references within a function
struct MmdbRef {
  unsigned line;
  std::string type;
};

// Structure to hold information about a function with mmdb references
struct FunctionInfo {
  std::string name;
  std::string file;
  unsigned startLine;
  unsigned endLine;
  std::set<std::string> mmdbTypes;
  std::vector<MmdbRef> mmdbRefs;
};

// Visitor to find mmdb namespace references
class MmdbFinderVisitor : public RecursiveASTVisitor<MmdbFinderVisitor> {
private:
  std::map<const FunctionDecl*, FunctionInfo> functionMap;
  const FunctionDecl* currentFunction = nullptr;
  ASTContext* context = nullptr;

public:
  MmdbFinderVisitor(ASTContext* ctx) : context(ctx) {}

  bool VisitFunctionDecl(FunctionDecl* func) {
    // Skip function templates and declarations without body
    if (func->isTemplated() || !func->hasBody()) {
      return true;
    }

    // Only process functions defined in the main source file, not in headers/libraries
    SourceManager& SM = context->getSourceManager();
    if (!SM.isInMainFile(func->getBeginLoc())) {
      return true;
    }

    currentFunction = func;
    FunctionInfo info;
    info.name = func->getNameAsString();
    info.file = SM.getFilename(func->getBeginLoc()).str();
    info.startLine = SM.getSpellingLineNumber(func->getBeginLoc());
    info.endLine = SM.getSpellingLineNumber(func->getEndLoc());

    // Store function info and visit its body
    functionMap[func] = info;
    TraverseStmt(func->getBody());
    currentFunction = nullptr;

    return true;
  }

  bool VisitMemberExpr(MemberExpr* expr) {
    if (currentFunction && functionMap.count(currentFunction)) {
      checkForMmdbNamespace(expr, functionMap[currentFunction]);
    }
    return true;
  }

  bool VisitDeclRefExpr(DeclRefExpr* expr) {
    if (currentFunction && functionMap.count(currentFunction)) {
      // Check if this is a reference to an mmdb type or function
      std::string name = expr->getNameInfo().getAsString();
      if (name.find("mmdb::") != std::string::npos) {
        unsigned line = context->getSourceManager().getSpellingLineNumber(expr->getBeginLoc());
        functionMap[currentFunction].mmdbTypes.insert(name);
        functionMap[currentFunction].mmdbRefs.push_back({line, name});
      }
    }
    return true;
  }

  bool VisitCallExpr(CallExpr* expr) {
    if (currentFunction && functionMap.count(currentFunction)) {
      // Get the called function/method
      if (FunctionDecl* calledFunc = expr->getDirectCallee()) {
        std::string funcName = calledFunc->getQualifiedNameAsString();
        if (funcName.find("mmdb::") != std::string::npos) {
          unsigned line = context->getSourceManager().getSpellingLineNumber(expr->getBeginLoc());
          functionMap[currentFunction].mmdbTypes.insert(funcName);
          functionMap[currentFunction].mmdbRefs.push_back({line, funcName});
        }
      }
    }
    return true;
  }

  bool VisitTypeLoc(TypeLoc tl) {
    if (currentFunction && functionMap.count(currentFunction)) {
      // Check type names for mmdb references
      std::string typeStr = tl.getType().getAsString();
      if (typeStr.find("mmdb::") != std::string::npos) {
        unsigned line = context->getSourceManager().getSpellingLineNumber(tl.getBeginLoc());
        functionMap[currentFunction].mmdbTypes.insert(typeStr);
        functionMap[currentFunction].mmdbRefs.push_back({line, typeStr});
      }
    }
    return true;
  }

private:
  void checkForMmdbNamespace(MemberExpr* expr, FunctionInfo& info) {
    // Get the full qualified name
    std::string name = expr->getMemberNameInfo().getAsString();

    // Check base expression for mmdb namespace
    Expr* base = expr->getBase()->IgnoreParenImpCasts();
    std::string baseType = base->getType().getAsString();

    if (baseType.find("mmdb::") != std::string::npos || name.find("mmdb::") != std::string::npos) {
      unsigned line = context->getSourceManager().getSpellingLineNumber(expr->getBeginLoc());
      info.mmdbTypes.insert(baseType);
      info.mmdbRefs.push_back({line, baseType});
    }
  }

public:
  const std::map<const FunctionDecl*, FunctionInfo>& getFunctionMap() const {
    return functionMap;
  }
};

// Consumer to process the AST
class MmdbConsumer : public ASTConsumer {
private:
  std::vector<FunctionInfo> functionsWithMmdb;
  ASTContext* context = nullptr;

public:
  virtual ~MmdbConsumer() = default;

  void HandleTranslationUnit(ASTContext& ctx) override {
    context = &ctx;
    MmdbFinderVisitor visitor(&ctx);
    visitor.TraverseDecl(ctx.getTranslationUnitDecl());

    // Filter to only functions that have mmdb references
    for (const auto& pair : visitor.getFunctionMap()) {
      const FunctionInfo& info = pair.second;
      if (!info.mmdbTypes.empty()) {
        functionsWithMmdb.push_back(info);
      }
    }
  }

  const std::vector<FunctionInfo>& getFunctionsWithMmdb() const {
    return functionsWithMmdb;
  }
};

// Global storage for results
std::vector<FunctionInfo> g_allFunctions;
std::set<std::string> g_allFiles;

// Map from source file path -> ordered list of #include spellings ("header.h", <header.h>, etc.)
// We use a map-of-set to deduplicate while still preserving insertion order in a
// companion vector. Simplest: set for dedup + vector for order.
struct IncludeList {
  std::vector<std::string> order;
  std::set<std::string> seen;
  void add(const std::string& inc) {
    if (seen.insert(inc).second) order.push_back(inc);
  }
};
std::map<std::string, IncludeList> g_includesByFile;

// PPCallbacks that records every #include directive seen in the main source
// file (not in nested headers), along with how it was spelled in the source
// so we can reproduce it exactly.
class IncludeCapture : public PPCallbacks {
public:
  IncludeCapture(SourceManager& SM) : SM_(SM) {}

  void InclusionDirective(SourceLocation HashLoc,
                          const Token& /*IncludeTok*/,
                          StringRef FileName,
                          bool IsAngled,
                          CharSourceRange /*FilenameRange*/,
                          OptionalFileEntryRef /*File*/,
                          StringRef /*SearchPath*/,
                          StringRef /*RelativePath*/,
                          const Module* /*SuggestedModule*/,
                          bool /*ModuleImported*/,
                          SrcMgr::CharacteristicKind /*FileType*/) override {
    // Only capture includes that appear directly in the translation unit's
    // main source file — skip anything transitively pulled in.
    if (!SM_.isInMainFile(HashLoc)) return;
    FileID mainID = SM_.getMainFileID();
    const FileEntry* mainFE = SM_.getFileEntryForID(mainID);
    if (!mainFE) return;
    std::string mainPath = mainFE->tryGetRealPathName().str();
    if (mainPath.empty()) {
      // tryGetRealPathName may be empty for in-memory files; fall back to
      // the filename the SourceManager knows about for the main file.
      mainPath = SM_.getFilename(SM_.getLocForStartOfFile(mainID)).str();
    }

    std::string spelled = IsAngled
      ? "<" + FileName.str() + ">"
      : "\"" + FileName.str() + "\"";
    g_includesByFile[mainPath].add(spelled);
  }

private:
  SourceManager& SM_;
};

// Custom consumer that stores results when destroyed
class MmdbConsumerWithCleanup : public MmdbConsumer {
public:
  ~MmdbConsumerWithCleanup() override {
    // Collect results when consumer is destroyed
    const auto& functions = getFunctionsWithMmdb();
    for (const auto& func : functions) {
      g_allFunctions.push_back(func);
      g_allFiles.insert(func.file);
    }
  }
};

// Action to create the consumer
class MmdbAction : public ASTFrontendAction {
public:
  std::unique_ptr<ASTConsumer> CreateASTConsumer(CompilerInstance& CI,
                                                  StringRef file) override {
    // Register a preprocessor callback to capture every #include directive
    // seen in the main source file of this translation unit.
    Preprocessor& pp = CI.getPreprocessor();
    pp.addPPCallbacks(std::make_unique<IncludeCapture>(CI.getSourceManager()));
    return std::make_unique<MmdbConsumerWithCleanup>();
  }
};

// Factory to create MmdbAction instances
class MmdbActionFactory : public FrontendActionFactory {
public:
  std::unique_ptr<FrontendAction> create() override {
    return std::make_unique<MmdbAction>();
  }
};

int main(int argc, const char** argv) {
  llvm::cl::OptionCategory tool_category("find-mmdb-functions options");
  llvm::cl::opt<std::string> coot_root(
      "coot-root", llvm::cl::desc("Root directory of coot source"),
      llvm::cl::init("/Users/dialpuri/lmb/coot"),
      llvm::cl::cat(tool_category));
  llvm::cl::opt<std::string> output_file(
      "output", llvm::cl::desc("Output JSON file"), llvm::cl::init("mmdb_functions.json"),
      llvm::cl::cat(tool_category));

  auto expected_parser = CommonOptionsParser::create(argc, argv, tool_category);
  if (!expected_parser) {
    llvm::errs() << expected_parser.takeError();
    return 1;
  }

  CommonOptionsParser& options_parser = expected_parser.get();
  ClangTool tool(options_parser.getCompilations(),
                 options_parser.getSourcePathList());

  // Inject clang's resource directory so built-in headers (stdarg.h etc.) are found
  std::string resource_dir;
  FILE* pipe = popen("clang -print-resource-dir 2>/dev/null", "r");
  if (pipe) {
    char buf[512];
    if (fgets(buf, sizeof(buf), pipe)) {
      resource_dir = buf;
      // Strip trailing newline
      if (!resource_dir.empty() && resource_dir.back() == '\n')
        resource_dir.pop_back();
    }
    pclose(pipe);
  }

  if (resource_dir.empty()) {
    llvm::errs() << "Warning: could not determine clang resource directory\n";
  } else {
    llvm::outs() << "Using clang resource dir: " << resource_dir << "\n";
    tool.appendArgumentsAdjuster(
        [resource_dir](const CommandLineArguments& args, StringRef /*Filename*/) {
          CommandLineArguments adjusted = args;
          adjusted.push_back("-resource-dir");
          adjusted.push_back(resource_dir);
          return adjusted;
        });
  }


  // Run the tool
  MmdbActionFactory factory;
  int result = tool.run(&factory);
  if (result != 0) {
    llvm::errs() << "Tool execution had errors\n";
  }

  const std::string COOT_ROOT = coot_root.getValue();

  // Group functions by file
  std::map<std::string, std::vector<const FunctionInfo*>> byFile;
  int totalRefs = 0;
  for (const auto& func : g_allFunctions) {
    byFile[func.file].push_back(&func);
    totalRefs += static_cast<int>(func.mmdbRefs.size());
  }

  // Build JSON output in the format expected by the visualiser
  llvm::json::Object output;

  // meta
  llvm::json::Object meta;
  meta["root"] = COOT_ROOT;
  meta["total_files_scanned"] = static_cast<int>(g_allFiles.size());
  meta["files_with_mmdb_refs"] = static_cast<int>(byFile.size());
  meta["total_mmdb_refs"] = totalRefs;
  output["meta"] = std::move(meta);

  // files[] — one entry per source file, with nested functions[]
  llvm::json::Array filesArray;
  for (const auto& [filePath, funcs] : byFile) {
    llvm::json::Object fileObj;

    // Compute rel_path by stripping the coot root prefix
    std::string relPath = filePath;
    if (relPath.rfind(COOT_ROOT, 0) == 0)
      relPath = relPath.substr(std::string(COOT_ROOT).size() + 1);

    fileObj["rel_path"] = relPath;
    fileObj["includes_mmdb"] = true;

    // All #include directives captured from this file's main translation unit,
    // in source order, each spelled as it appears in the source
    // (e.g. "\"coot-utils/coot-coord-utils.hh\"" or "<mmdb2/mmdb_manager.h>").
    llvm::json::Array includesArray;
    llvm::json::Array mmdbIncludesArray;
    auto it = g_includesByFile.find(filePath);
    if (it != g_includesByFile.end()) {
      for (const auto& inc : it->second.order) {
        includesArray.push_back(inc);
        // Convenience subset: anything that mentions mmdb2 / mmdb_.
        if (inc.find("mmdb") != std::string::npos)
          mmdbIncludesArray.push_back(inc);
      }
    }
    fileObj["includes"]       = std::move(includesArray);
    fileObj["mmdb_includes"]  = std::move(mmdbIncludesArray);

    int fileTotalRefs = 0;
    llvm::json::Array functionsArray;
    for (const FunctionInfo* func : funcs) {
      llvm::json::Object funcObj;
      funcObj["name"] = func->name;
      funcObj["line"] = static_cast<int>(func->startLine);
      funcObj["end_line"] = static_cast<int>(func->endLine);
      funcObj["mmdb_ref_count"] = static_cast<int>(func->mmdbRefs.size());
      fileTotalRefs += static_cast<int>(func->mmdbRefs.size());

      llvm::json::Array symbolsArray;
      for (const auto& sym : func->mmdbTypes)
        symbolsArray.push_back(sym);
      funcObj["mmdb_symbols"] = std::move(symbolsArray);

      // refs as array of [line, type] pairs to match FunctionRecord
      llvm::json::Array refsArray;
      for (const auto& ref : func->mmdbRefs) {
        llvm::json::Array pair;
        pair.push_back(static_cast<int>(ref.line));
        pair.push_back(ref.type);
        refsArray.push_back(std::move(pair));
      }
      funcObj["refs"] = std::move(refsArray);

      functionsArray.push_back(std::move(funcObj));
    }

    fileObj["functions"] = std::move(functionsArray);
    fileObj["total_mmdb_refs"] = fileTotalRefs;
    fileObj["file_level_refs"] = llvm::json::Array{};

    filesArray.push_back(std::move(fileObj));
  }
  output["files"] = std::move(filesArray);

  // Write output
  std::ofstream outfile(output_file);
  if (!outfile) {
    llvm::errs() << "Could not open output file: " << output_file << "\n";
    return 1;
  }

  std::string json_buf;
  llvm::raw_string_ostream json_stream(json_buf);
  json_stream << llvm::json::Value(std::move(output));
  json_stream.flush();

  outfile << json_buf << "\n";
  outfile.close();

  llvm::outs() << "Results written to " << output_file << "\n";
  llvm::outs() << "Found " << g_allFunctions.size() << " functions with mmdb references\n";
  llvm::outs() << "Across " << g_allFiles.size() << " files\n";

  return 0;
}
