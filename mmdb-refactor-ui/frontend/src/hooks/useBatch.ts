import { useState, useRef, useCallback, useEffect } from 'react'
import { generateAllTests } from '../api'
import type { BatchEvent } from '../api'

export interface BatchLogEntry {
  type: 'done' | 'skip' | 'error'
  fn: string
  file: string
  message?: string
  // done-specific
  mmdb_status?: string
  mmdb_attempts?: number
  gemmi_status?: string
  gemmi_attempts?: number
}

export type BatchPhase = 'idle' | 'running' | 'finished'

export interface BatchState {
  phase: BatchPhase
  total: number
  done: number
  skipped: number
  errors: number
  currentFn: string | null
  currentFile: string | null
  currentAttempt: { variant: string; attempt: number; max: number; error: string } | null
  log: BatchLogEntry[]
}

export interface BatchConfig {
  model: string
  skipExisting: boolean
  extraInstructions: string
}

const IDLE: BatchState = {
  phase: 'idle', total: 0, done: 0, skipped: 0, errors: 0,
  currentFn: null, currentFile: null, currentAttempt: null, log: [],
}

export function useBatch(onFinish: () => void) {
  const onFinishRef = useRef(onFinish)
  useEffect(() => { onFinishRef.current = onFinish }, [onFinish])

  const [config, setConfig] = useState<BatchConfig>({
    model: 'gemma4',
    skipExisting: true,
    extraInstructions: '',
  })
  const [state, setState] = useState<BatchState>(IDLE)
  const abortRef = useRef<AbortController | null>(null)

  const start = useCallback(async () => {
    setState({ ...IDLE, phase: 'running' })
    abortRef.current = new AbortController()

    try {
      const gen = generateAllTests(
        config.model,
        config.skipExisting,
        config.extraInstructions,
        abortRef.current.signal,
      )

      for await (const ev of gen) {
        setState(prev => applyEvent(prev, ev))
      }
    } catch (e: unknown) {
      if (!(e instanceof Error && e.name === 'AbortError')) {
        setState(prev => ({
          ...prev,
          log: [...prev.log, { type: 'error', fn: '', file: '', message: String(e) }],
        }))
      }
    } finally {
      setState(prev => ({ ...prev, phase: 'finished', currentFn: null, currentFile: null, currentAttempt: null }))
      onFinishRef.current()
    }
  }, [config])

  const stop = useCallback(() => {
    abortRef.current?.abort()
  }, [])

  const reset = useCallback(() => {
    setState(IDLE)
  }, [])

  return { config, setConfig, state, start, stop, reset }
}

function applyEvent(prev: BatchState, ev: BatchEvent): BatchState {
  switch (ev.type) {
    case 'start':
      return { ...prev, total: ev.total }
    case 'progress':
      return { ...prev, done: ev.done, skipped: ev.skipped, errors: ev.errors, currentFn: ev.fn, currentFile: ev.file, currentAttempt: null }
    case 'skip':
      return { ...prev, done: ev.done, skipped: ev.skipped, errors: ev.errors, currentAttempt: null,
               log: [...prev.log, { type: 'skip', fn: ev.fn, file: ev.file }] }
    case 'attempt':
      return { ...prev, done: ev.done, skipped: ev.skipped, errors: ev.errors,
               currentAttempt: { variant: ev.variant, attempt: ev.attempt, max: ev.max, error: ev.error } }
    case 'done':
      return { ...prev, done: ev.done, skipped: ev.skipped, errors: ev.errors, currentAttempt: null,
               log: [...prev.log, { type: 'done', fn: ev.fn, file: ev.file,
                 mmdb_status: ev.mmdb_status, mmdb_attempts: ev.mmdb_attempts,
                 gemmi_status: ev.gemmi_status, gemmi_attempts: ev.gemmi_attempts }] }
    case 'error':
      return { ...prev, done: ev.done, skipped: ev.skipped, errors: ev.errors, currentAttempt: null,
               log: [...prev.log, { type: 'error', fn: ev.fn, file: ev.file, message: ev.message }] }
    case 'finish':
      return { ...prev, done: ev.done, skipped: ev.skipped, errors: ev.errors, total: ev.total }
    default:
      return prev
  }
}
