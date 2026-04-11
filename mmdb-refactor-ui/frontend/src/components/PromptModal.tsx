import { useEffect, useState } from 'react'

export interface ModalSection {
  title: string
  text: string
  /** If true, this section is expanded by default. Defaults to last section. */
  defaultOpen?: boolean
}

interface Props {
  title: string
  sections: ModalSection[] | null
  loading: boolean
  error: string | null
  onClose: () => void
}

function CopyButton({ text }: { text: string }) {
  const [copied, setCopied] = useState(false)
  const copy = () => {
    navigator.clipboard.writeText(text).then(() => {
      setCopied(true)
      setTimeout(() => setCopied(false), 1500)
    })
  }
  return (
    <button
      onClick={copy}
      className="text-xs text-zinc-500 hover:text-zinc-300 transition-colors px-2 py-0.5 rounded border border-zinc-700 hover:border-zinc-500"
    >
      {copied ? 'Copied!' : 'Copy'}
    </button>
  )
}

function Section({ title, text, defaultOpen = false }: ModalSection) {
  const [open, setOpen] = useState(defaultOpen)
  return (
    <div className="border border-zinc-800 rounded">
      <button
        className="w-full flex items-center justify-between px-3 py-2 text-xs font-semibold text-zinc-400 hover:text-zinc-200 bg-zinc-900 rounded-t transition-colors"
        onClick={() => setOpen(o => !o)}
      >
        <span>{title}</span>
        <div className="flex items-center gap-2">
          {open && <CopyButton text={text} />}
          <span className="text-zinc-600">{open ? '▲' : '▼'}</span>
        </div>
      </button>
      {open && (
        <pre className="p-3 text-xs font-mono text-zinc-300 bg-zinc-950 whitespace-pre-wrap overflow-auto max-h-80 rounded-b leading-5">
          {text || <span className="text-zinc-600 italic">empty</span>}
        </pre>
      )}
    </div>
  )
}

export default function PromptModal({ title, sections, loading, error, onClose }: Props) {
  useEffect(() => {
    const handler = (e: KeyboardEvent) => { if (e.key === 'Escape') onClose() }
    window.addEventListener('keydown', handler)
    return () => window.removeEventListener('keydown', handler)
  }, [onClose])

  return (
    <div
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/70"
      onMouseDown={e => { if (e.target === e.currentTarget) onClose() }}
    >
      <div className="bg-zinc-900 border border-zinc-700 rounded-lg shadow-2xl w-full max-w-3xl mx-4 flex flex-col max-h-[85vh]">
        {/* Header */}
        <div className="flex items-center justify-between px-4 py-3 border-b border-zinc-800 flex-shrink-0">
          <span className="text-sm font-semibold text-zinc-200">{title}</span>
          <button
            onClick={onClose}
            className="text-zinc-500 hover:text-zinc-200 text-lg leading-none transition-colors"
          >
            ×
          </button>
        </div>

        {/* Body */}
        <div className="flex-1 overflow-auto p-4 space-y-3">
          {loading && (
            <div className="text-xs text-zinc-500 animate-pulse">Loading…</div>
          )}
          {error && (
            <div className="text-xs text-red-400 bg-red-950/40 border border-red-900/50 rounded px-3 py-2">{error}</div>
          )}
          {sections && sections.map((s, i) => (
            <Section
              key={s.title}
              title={s.title}
              text={s.text}
              defaultOpen={s.defaultOpen ?? i === sections.length - 1}
            />
          ))}
        </div>

        {/* Footer */}
        <div className="px-4 py-3 border-t border-zinc-800 flex justify-end flex-shrink-0">
          <button onClick={onClose} className="btn btn-secondary text-sm">Close</button>
        </div>
      </div>
    </div>
  )
}
