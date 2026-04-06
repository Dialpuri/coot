import { useState, useRef, useCallback, useEffect } from 'react'

/**
 * Drag a bottom handle DOWN to grow, UP to shrink.
 * Returns [size, onMouseDown].
 */
export function useVerticalResize(initial: number, min: number, max: number) {
  const [size, setSize] = useState(initial)

  const startDrag = useCallback(
    (e: React.MouseEvent) => {
      e.preventDefault()
      const startY = e.clientY
      const startSize = size

      const onMove = (ev: MouseEvent) => {
        const delta = ev.clientY - startY
        setSize(Math.max(min, Math.min(max, startSize + delta)))
      }
      const onUp = () => {
        window.removeEventListener('mousemove', onMove)
        window.removeEventListener('mouseup', onUp)
      }
      window.addEventListener('mousemove', onMove)
      window.addEventListener('mouseup', onUp)
    },
    [size, min, max],
  )

  return [size, setSize, startDrag] as const
}

/**
 * Drag a top handle UP to grow, DOWN to shrink (inverted direction).
 * Returns [size, onMouseDown].
 */
export function useVerticalResizeInverted(initial: number, min: number, max: number) {
  const [size, setSize] = useState(initial)

  const startDrag = useCallback(
    (e: React.MouseEvent) => {
      e.preventDefault()
      const startY = e.clientY
      const startSize = size

      const onMove = (ev: MouseEvent) => {
        const delta = startY - ev.clientY  // inverted: up = bigger
        setSize(Math.max(min, Math.min(max, startSize + delta)))
      }
      const onUp = () => {
        window.removeEventListener('mousemove', onMove)
        window.removeEventListener('mouseup', onUp)
      }
      window.addEventListener('mousemove', onMove)
      window.addEventListener('mouseup', onUp)
    },
    [size, min, max],
  )

  return [size, setSize, startDrag] as const
}

/**
 * Drag a right-edge handle RIGHT to grow, LEFT to shrink.
 * Optionally persists the size to localStorage via `storageKey`.
 * Returns [size, onMouseDown].
 */
export function useHorizontalResize(
  initial: number,
  min: number,
  max: number,
  storageKey?: string,
) {
  const [size, setSize] = useState(() => {
    if (storageKey) {
      const saved = localStorage.getItem(storageKey)
      if (saved !== null) {
        const n = parseInt(saved, 10)
        if (!isNaN(n) && n >= min && n <= max) return n
      }
    }
    return initial
  })

  useEffect(() => {
    if (storageKey) localStorage.setItem(storageKey, String(size))
  }, [size, storageKey])

  const startDrag = useCallback(
    (e: React.MouseEvent) => {
      e.preventDefault()
      const startX = e.clientX
      const startSize = size

      const onMove = (ev: MouseEvent) => {
        const delta = ev.clientX - startX
        setSize(Math.max(min, Math.min(max, startSize + delta)))
      }
      const onUp = () => {
        window.removeEventListener('mousemove', onMove)
        window.removeEventListener('mouseup', onUp)
      }
      window.addEventListener('mousemove', onMove)
      window.addEventListener('mouseup', onUp)
    },
    [size, min, max],
  )

  return [size, startDrag] as const
}

/**
 * Drag a vertical divider to adjust a left/right percentage split.
 * Returns [leftPercent, containerRef, onMouseDown].
 */
export function useHorizontalSplit(initial = 50) {
  const [leftPercent, setLeftPercent] = useState(initial)
  const containerRef = useRef<HTMLDivElement>(null)

  const startDrag = useCallback((e: React.MouseEvent) => {
    e.preventDefault()

    const onMove = (ev: MouseEvent) => {
      if (!containerRef.current) return
      const rect = containerRef.current.getBoundingClientRect()
      const pct = ((ev.clientX - rect.left) / rect.width) * 100
      setLeftPercent(Math.max(15, Math.min(85, pct)))
    }
    const onUp = () => {
      window.removeEventListener('mousemove', onMove)
      window.removeEventListener('mouseup', onUp)
    }
    window.addEventListener('mousemove', onMove)
    window.addEventListener('mouseup', onUp)
  }, [])

  return [leftPercent, containerRef, startDrag] as const
}
