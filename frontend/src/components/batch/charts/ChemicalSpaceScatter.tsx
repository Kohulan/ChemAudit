/**
 * ChemicalSpaceScatter (VIZ-06)
 *
 * Canvas 2D chemical space scatter plot for PCA/t-SNE coordinates.
 * Supports >800 points without browser jank.
 * Features: color-by property, hit detection, brush selection rectangle,
 * click-to-toggle, and PNG download via canvas.toBlob().
 */

import React, { useRef, useEffect, useCallback, useState, useMemo } from 'react';
import { batchApi } from '../../../services/api';
import type { ChemSpaceCoordinates } from '../../../types/analytics';
import type { BatchResult } from '../../../types/batch';

interface ChemicalSpaceScatterProps {
  data: ChemSpaceCoordinates | null;
  results: BatchResult[];
  colorByProperty: string;
  onColorByChange: (prop: string) => void;
  selectedIndices: Set<number>;
  onSelectionChange: (indices: Set<number>) => void;
  jobId: string;
}

const PROPERTY_OPTIONS = ['overall_score', 'QED', 'SA_score', 'Fsp3'];
const MARGIN = 40;
const CANVAS_HEIGHT = 400;

function getPropertyValue(r: BatchResult, prop: string): number | null {
  if (r.status !== 'success') return null;
  switch (prop) {
    case 'overall_score':
      return r.validation?.overall_score ?? null;
    case 'QED':
      return r.scoring?.druglikeness?.qed_score ?? null;
    case 'SA_score':
      return r.scoring?.admet?.sa_score ?? null;
    case 'Fsp3':
      return r.scoring?.admet?.fsp3 ?? null;
    default:
      return null;
  }
}

function valueToColor(value: number, min: number, max: number): string {
  const t = max > min ? Math.max(0, Math.min(1, (value - min) / (max - min))) : 0.5;
  // Blue (low) -> Green (mid) -> Yellow (high)
  const r = Math.round(255 * t);
  const g = Math.round(200 * (1 - Math.abs(t - 0.5) * 2) + 55);
  const b = Math.round(255 * (1 - t));
  return `rgb(${r}, ${g}, ${b})`;
}

interface PointData {
  px: number;
  py: number;
  index: number;
  color: string;
}

export const ChemicalSpaceScatter = React.memo(function ChemicalSpaceScatter({
  data,
  results,
  colorByProperty,
  onColorByChange,
  selectedIndices,
  onSelectionChange,
  jobId,
}: ChemicalSpaceScatterProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);
  const pointsRef = useRef<PointData[]>([]);
  const [tooltipState, setTooltipState] = useState<{
    x: number;
    y: number;
    index: number;
    smiles: string;
  } | null>(null);
  const [canvasWidth, setCanvasWidth] = useState(600);
  const [brushRect, setBrushRect] = useState<{
    startX: number;
    startY: number;
    endX: number;
    endY: number;
  } | null>(null);
  const brushStartRef = useRef<{ x: number; y: number } | null>(null);
  const [triggeringTsne, setTriggeringTsne] = useState(false);

  // Build color value range
  const { colorMin, colorMax, colorMap } = useMemo(() => {
    let cMin = Infinity;
    let cMax = -Infinity;
    const map = new Map<number, number>();
    for (const r of results) {
      const val = getPropertyValue(r, colorByProperty);
      if (val !== null) {
        map.set(r.index, val);
        if (val < cMin) cMin = val;
        if (val > cMax) cMax = val;
      }
    }
    return {
      colorMin: cMin === Infinity ? 0 : cMin,
      colorMax: cMax === -Infinity ? 1 : cMax,
      colorMap: map,
    };
  }, [results, colorByProperty]);

  // ResizeObserver for responsive width
  useEffect(() => {
    const container = containerRef.current;
    if (!container) return;

    const observer = new ResizeObserver((entries) => {
      for (const entry of entries) {
        setCanvasWidth(entry.contentRect.width);
      }
    });
    observer.observe(container);
    return () => observer.disconnect();
  }, []);

  // Draw canvas
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !data) return;

    const dpr = window.devicePixelRatio || 1;
    const logicalW = canvasWidth;
    const logicalH = CANVAS_HEIGHT;

    canvas.width = logicalW * dpr;
    canvas.height = logicalH * dpr;
    canvas.style.width = `${logicalW}px`;
    canvas.style.height = `${logicalH}px`;

    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    ctx.scale(dpr, dpr);
    ctx.clearRect(0, 0, logicalW, logicalH);

    // Background
    ctx.fillStyle = getComputedStyle(document.documentElement)
      .getPropertyValue('--color-surface-sunken')
      .trim() || '#f8f9fa';
    ctx.fillRect(0, 0, logicalW, logicalH);

    if (data.coordinates.length === 0) return;

    // Compute coordinate range
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    for (const coord of data.coordinates) {
      if (coord[0] < xMin) xMin = coord[0];
      if (coord[0] > xMax) xMax = coord[0];
      if (coord[1] < yMin) yMin = coord[1];
      if (coord[1] > yMax) yMax = coord[1];
    }
    const xRange = xMax - xMin || 1;
    const yRange = yMax - yMin || 1;
    const plotW = logicalW - 2 * MARGIN;
    const plotH = logicalH - 2 * MARGIN;

    // Draw grid lines
    ctx.strokeStyle = 'rgba(128, 128, 128, 0.15)';
    ctx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const x = MARGIN + (plotW * i) / 4;
      const y = MARGIN + (plotH * i) / 4;
      ctx.beginPath();
      ctx.moveTo(x, MARGIN);
      ctx.lineTo(x, MARGIN + plotH);
      ctx.stroke();
      ctx.beginPath();
      ctx.moveTo(MARGIN, y);
      ctx.lineTo(MARGIN + plotW, y);
      ctx.stroke();
    }

    // Draw axes labels
    ctx.fillStyle = 'rgba(128, 128, 128, 0.8)';
    ctx.font = '11px sans-serif';
    ctx.textAlign = 'center';
    const xLabel = data.method === 'pca' ? 'PC1' : 't-SNE 1';
    const yLabel = data.method === 'pca' ? 'PC2' : 't-SNE 2';
    ctx.fillText(xLabel, logicalW / 2, logicalH - 8);
    ctx.save();
    ctx.translate(12, logicalH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(yLabel, 0, 0);
    ctx.restore();

    // Map and draw points
    const hasSelection = selectedIndices.size > 0;
    const newPoints: PointData[] = [];

    for (let i = 0; i < data.coordinates.length; i++) {
      const [cx, cy] = data.coordinates[i];
      const molIdx = data.molecule_indices[i];
      const px = MARGIN + ((cx - xMin) / xRange) * plotW;
      const py = MARGIN + plotH - ((cy - yMin) / yRange) * plotH;

      const cVal = colorMap.get(molIdx) ?? 0;
      const color = valueToColor(cVal, colorMin, colorMax);

      newPoints.push({ px, py, index: molIdx, color });
    }

    // Draw non-selected points first, then selected on top
    for (const pt of newPoints) {
      const isSelected = selectedIndices.has(pt.index);
      if (isSelected) continue;

      ctx.globalAlpha = hasSelection ? 0.3 : 0.7;
      ctx.beginPath();
      ctx.arc(pt.px, pt.py, 4, 0, Math.PI * 2);
      ctx.fillStyle = pt.color;
      ctx.fill();
    }

    for (const pt of newPoints) {
      if (!selectedIndices.has(pt.index)) continue;

      ctx.globalAlpha = 1;
      ctx.beginPath();
      ctx.arc(pt.px, pt.py, 6, 0, Math.PI * 2);
      ctx.fillStyle = '#c41e3a';
      ctx.fill();
      ctx.strokeStyle = '#fff';
      ctx.lineWidth = 1.5;
      ctx.stroke();
    }

    ctx.globalAlpha = 1;

    // Draw brush selection rectangle if active
    if (brushRect) {
      ctx.strokeStyle = 'rgba(37, 99, 235, 0.7)';
      ctx.lineWidth = 1.5;
      ctx.setLineDash([6, 4]);
      ctx.strokeRect(
        brushRect.startX,
        brushRect.startY,
        brushRect.endX - brushRect.startX,
        brushRect.endY - brushRect.startY
      );
      ctx.fillStyle = 'rgba(37, 99, 235, 0.08)';
      ctx.fillRect(
        brushRect.startX,
        brushRect.startY,
        brushRect.endX - brushRect.startX,
        brushRect.endY - brushRect.startY
      );
      ctx.setLineDash([]);
    }

    pointsRef.current = newPoints;
  }, [data, canvasWidth, selectedIndices, colorMap, colorMin, colorMax, brushRect]);

  // Hit detection on mousemove (debounced)
  const mouseMoveTimer = useRef<ReturnType<typeof setTimeout> | null>(null);

  const handleMouseMove = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      const canvas = canvasRef.current;
      if (!canvas) return;

      // If brushing, update the brush rectangle
      if (brushStartRef.current) {
        const rect = canvas.getBoundingClientRect();
        const x = e.clientX - rect.left;
        const y = e.clientY - rect.top;
        setBrushRect({
          startX: brushStartRef.current.x,
          startY: brushStartRef.current.y,
          endX: x,
          endY: y,
        });
        return;
      }

      if (mouseMoveTimer.current) clearTimeout(mouseMoveTimer.current);

      mouseMoveTimer.current = setTimeout(() => {
        const rect = canvas.getBoundingClientRect();
        const mx = e.clientX - rect.left;
        const my = e.clientY - rect.top;

        let nearest: PointData | null = null;
        let nearestDist = 12;

        for (const pt of pointsRef.current) {
          const dx = pt.px - mx;
          const dy = pt.py - my;
          const dist = Math.sqrt(dx * dx + dy * dy);
          if (dist < nearestDist) {
            nearest = pt;
            nearestDist = dist;
          }
        }

        if (nearest) {
          const r = results.find((r) => r.index === nearest!.index);
          setTooltipState({
            x: e.clientX,
            y: e.clientY,
            index: nearest.index,
            smiles: r?.smiles || `Molecule ${nearest.index}`,
          });
        } else {
          setTooltipState(null);
        }
      }, 50);
    },
    [results]
  );

  const handleMouseLeave = useCallback(() => {
    setTooltipState(null);
    if (mouseMoveTimer.current) clearTimeout(mouseMoveTimer.current);
  }, []);

  // Click to toggle point
  const handleClick = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      // Don't toggle if we just finished a brush
      if (brushStartRef.current) return;

      const canvas = canvasRef.current;
      if (!canvas) return;

      const rect = canvas.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      let nearest: PointData | null = null;
      let nearestDist = 12;

      for (const pt of pointsRef.current) {
        const dx = pt.px - mx;
        const dy = pt.py - my;
        const dist = Math.sqrt(dx * dx + dy * dy);
        if (dist < nearestDist) {
          nearest = pt;
          nearestDist = dist;
        }
      }

      if (nearest) {
        const next = new Set(selectedIndices);
        if (next.has(nearest.index)) {
          next.delete(nearest.index);
        } else {
          next.add(nearest.index);
        }
        onSelectionChange(next);
      }
    },
    [selectedIndices, onSelectionChange]
  );

  // Brush selection: mousedown starts the rectangle
  const handleMouseDown = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      const canvas = canvasRef.current;
      if (!canvas) return;

      const rect = canvas.getBoundingClientRect();
      brushStartRef.current = {
        x: e.clientX - rect.left,
        y: e.clientY - rect.top,
      };
    },
    []
  );

  // Brush selection: mouseup selects all points inside the rectangle
  const handleMouseUp = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      if (!brushStartRef.current) return;

      const canvas = canvasRef.current;
      if (!canvas) {
        brushStartRef.current = null;
        setBrushRect(null);
        return;
      }

      const rect = canvas.getBoundingClientRect();
      const endX = e.clientX - rect.left;
      const endY = e.clientY - rect.top;

      const sx = Math.min(brushStartRef.current.x, endX);
      const sy = Math.min(brushStartRef.current.y, endY);
      const ex = Math.max(brushStartRef.current.x, endX);
      const ey = Math.max(brushStartRef.current.y, endY);

      // Minimum brush size to distinguish from click
      const brushW = ex - sx;
      const brushH = ey - sy;

      if (brushW > 5 && brushH > 5) {
        const newSelection = new Set<number>();
        for (const pt of pointsRef.current) {
          if (pt.px >= sx && pt.px <= ex && pt.py >= sy && pt.py <= ey) {
            newSelection.add(pt.index);
          }
        }
        if (newSelection.size > 0) {
          onSelectionChange(newSelection);
        }
      }

      brushStartRef.current = null;
      setBrushRect(null);
    },
    [onSelectionChange]
  );

  // Compute t-SNE trigger
  const handleTriggerTsne = useCallback(async () => {
    setTriggeringTsne(true);
    try {
      await batchApi.triggerAnalytics(jobId, 'chemical_space', { method: 'tsne' });
    } catch {
      // Ignore — parent polling will pick it up
    }
    setTriggeringTsne(false);
  }, [jobId]);

  // PNG download
  const handleDownload = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    canvas.toBlob((blob) => {
      if (!blob) return;
      const a = document.createElement('a');
      a.href = URL.createObjectURL(blob);
      a.download = 'chemical-space.png';
      a.click();
      URL.revokeObjectURL(a.href);
    }, 'image/png');
  }, []);

  if (!data) {
    return (
      <div className="flex items-center justify-center h-[400px] text-sm text-[var(--color-text-muted)]">
        <div className="animate-pulse space-y-3 w-full">
          <div className="h-4 w-48 bg-[var(--color-border)] rounded mx-auto" />
          <div className="h-[360px] bg-[var(--color-border)]/50 rounded-xl" />
        </div>
      </div>
    );
  }

  const canShowTsne =
    data.method === 'pca' && results.length <= 2000;

  return (
    <div className="space-y-3">
      {/* Controls */}
      <div className="flex flex-wrap items-center gap-3">
        <label className="flex items-center gap-1.5 text-xs text-[var(--color-text-secondary)]">
          Color:
          <select
            value={colorByProperty}
            onChange={(e) => onColorByChange(e.target.value)}
            className="text-xs px-2 py-1 rounded-md bg-[var(--color-surface-sunken)] border border-[var(--color-border)] text-[var(--color-text-primary)]"
          >
            {PROPERTY_OPTIONS.map((p) => (
              <option key={p} value={p}>
                {p}
              </option>
            ))}
          </select>
        </label>

        {canShowTsne && (
          <button
            onClick={handleTriggerTsne}
            disabled={triggeringTsne}
            className="text-xs px-3 py-1 rounded-md bg-[var(--color-primary)]/10 text-[var(--color-primary)] hover:bg-[var(--color-primary)]/20 transition-colors disabled:opacity-50"
          >
            {triggeringTsne ? 'Computing...' : 'Compute t-SNE'}
          </button>
        )}

        <button
          onClick={handleDownload}
          className="text-xs px-3 py-1 rounded-md bg-[var(--color-surface-sunken)] text-[var(--color-text-secondary)] hover:text-[var(--color-text-primary)] transition-colors ml-auto"
        >
          Download PNG
        </button>

        <span className="text-xs text-[var(--color-text-muted)]">
          {data.method.toUpperCase()} | {data.coordinates.length} points
        </span>
      </div>

      {/* Canvas */}
      <div ref={containerRef} className="relative w-full rounded-xl overflow-hidden">
        <canvas
          ref={canvasRef}
          className="cursor-crosshair"
          onMouseMove={handleMouseMove}
          onMouseLeave={handleMouseLeave}
          onClick={handleClick}
          onMouseDown={handleMouseDown}
          onMouseUp={handleMouseUp}
        />

        {/* Tooltip */}
        {tooltipState && (
          <div
            className="fixed z-50 pointer-events-none rounded-lg bg-[var(--color-surface-elevated)] border border-[var(--color-border)] shadow-lg p-2 text-xs max-w-[300px]"
            style={{
              left: tooltipState.x + 12,
              top: tooltipState.y - 40,
            }}
          >
            <p className="font-mono text-[var(--color-text-primary)] truncate">
              {tooltipState.smiles.length > 40
                ? tooltipState.smiles.slice(0, 40) + '...'
                : tooltipState.smiles}
            </p>
            <p className="text-[var(--color-text-muted)]">Index: {tooltipState.index}</p>
          </div>
        )}
      </div>
    </div>
  );
});
