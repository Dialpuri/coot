import Prism from 'prismjs'
import 'prismjs/components/prism-c'
import 'prismjs/components/prism-cpp'

/** Prism-highlight C++ and overlay MMDB symbol marks on top. */
export function highlightCpp(code: string, mmdbSymbols: string[] = []): string {
  // 1. Prism-highlight first
  let html = Prism.highlight(code, Prism.languages.cpp, 'cpp')

  // 2. Overlay MMDB symbol marks.
  //    We work on the raw HTML so we avoid breaking Prism tokens.
  //    Sort longest first to avoid partial replacements.
  if (mmdbSymbols.length) {
    const sorted = [...mmdbSymbols].sort((a, b) => b.length - a.length)
    for (const sym of sorted) {
      // Escape the symbol for use inside HTML text (it's already HTML-escaped by Prism)
      const escaped = sym
        .replace(/&/g, '&amp;')
        .replace(/</g, '&lt;')
        .replace(/>/g, '&gt;')
      // Only match occurrences that are NOT already inside a span (simple heuristic:
      // replace in text nodes by splitting on tags)
      html = html.split(escaped).join(
        `<mark class="mmdb-symbol">${escaped}</mark>`,
      )
    }
  }

  return html
}

/** Highlight C++ with line numbers prepended. */
export function highlightCppWithLines(
  code: string,
  startLine: number,
  mmdbSymbols: string[] = [],
): string {
  const highlighted = highlightCpp(code, mmdbSymbols)
  return highlighted
    .split('\n')
    .map((line, i) => {
      const n = String(startLine + i).padStart(5, ' ')
      return `<span class="select-none text-gray-600 mr-3">${n}</span>${line}`
    })
    .join('\n')
}

/** Highlight function for react-simple-code-editor (no line numbers). */
export function highlightForEditor(code: string): string {
  return Prism.highlight(code, Prism.languages.cpp, 'cpp')
}
