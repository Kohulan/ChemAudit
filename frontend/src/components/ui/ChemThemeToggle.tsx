/**
 * ChemThemeToggle
 *
 * Chemistry-themed toggle switch for light/dark mode.
 * A round-bottom flask with liquid and bubbles slides between
 * a day scene (sun, clouds) and a night scene (moon, stars).
 */

import styled from 'styled-components';
import { useThemeContext } from '../../contexts/ThemeContext';
import { cn } from '../../lib/utils';

interface ChemThemeToggleProps {
  className?: string;
}

export function ChemThemeToggle({ className }: ChemThemeToggleProps) {
  const { resolvedTheme, setTheme } = useThemeContext();
  const isDark = resolvedTheme === 'dark';

  const handleChange = () => {
    setTheme(isDark ? 'light' : 'dark');
  };

  return (
    <StyledWrapper
      className={cn('overflow-visible transition-all duration-300', className)}
    >
      <label className="chem-toggle" aria-label={`Switch to ${isDark ? 'light' : 'dark'} mode`}>
        <input
          className="chem-toggle__checkbox"
          type="checkbox"
          checked={isDark}
          onChange={handleChange}
        />
        <div className="chem-toggle__container">
          <div className="chem-toggle__scenery">
            <div className="chem-toggle__star" />
            <div className="chem-toggle__star" />
            <div className="chem-toggle__star" />
            <div className="chem-toggle__star" />
            <div className="chem-toggle__star" />
            <div className="chem-toggle__star" />
            <div className="chem-toggle__star" />
            <div className="sun-primary" />
            <div className="sun-secondary" />
            <div className="moon" />
            <div className="moon-crater-1" />
            <div className="moon-crater-2" />
            <div className="chem-toggle__cloud" />
            <div className="chem-toggle__cloud" />
            <div className="chem-toggle__cloud" />
          </div>
          <div className="flask">
            <div className="flask__neck-container">
              <div className="flask__vapor" />
              <div className="flask__vapor" />
              <div className="flask__neck" />
            </div>
            <div className="flask__body" />
          </div>
          <div className="artificial__hidden">
            <div className="flask__shadow" />
          </div>
        </div>
      </label>
    </StyledWrapper>
  );
}

const StyledWrapper = styled.div`
  .chem-toggle {
    --toggle-size: 5.25px;
    --toggle-width: 10.625em;
    --toggle-height: 5.625em;
    --toggle-offset: calc((var(--toggle-height) - var(--flask-diameter)) / 2);
    --toggle-bg: linear-gradient(#2c4770, #070e2b 35%, #628cac 50% 70%, #a6c5d4)
      no-repeat;
    --flask-diameter: 4.375em;
    --radius: 99em;
    --transition: 0.4s;
    --liquid-day: #0d9488;
    --liquid-day-surface: rgba(13, 148, 136, 0.45);
    --liquid-night: #1e293b;
    --liquid-night-deep: #0f172a;
    --liquid-night-surface: rgba(30, 41, 59, 0.5);
    --glass: rgba(255, 255, 255, 0.82);
    --glass-edge: rgba(180, 210, 235, 0.45);
  }

  .chem-toggle,
  .chem-toggle *,
  .chem-toggle *::before,
  .chem-toggle *::after {
    box-sizing: border-box;
  }

  .chem-toggle {
    cursor: pointer;
    font-size: var(--toggle-size);
  }

  .chem-toggle__checkbox {
    display: none;
  }

  .chem-toggle__container {
    width: var(--toggle-width);
    height: var(--toggle-height);
    background: var(--toggle-bg);
    background-size: 100% 11.25em;
    background-position-y: -5.625em;
    border-radius: var(--radius);
    position: relative;
    transition: var(--transition);
    box-shadow:
      inset 1px 1px 3px rgba(255, 255, 255, 0.15),
      inset -1px -1px 3px rgba(0, 0, 0, 0.2),
      0 3px 10px rgba(0, 0, 0, 0.12),
      0 1px 3px rgba(0, 0, 0, 0.08);
    border: 2px solid rgba(255, 255, 255, 0.25);
  }

  /* ── Flask ── */

  .flask {
    display: flex;
    flex-direction: column;
    align-items: center;
    position: absolute;
    top: calc(var(--toggle-offset) - 1.688em + 0.188em);
    left: var(--toggle-offset);
    transition: var(--transition);
    z-index: 2;
  }

  .flask__neck-container {
    position: relative;
    z-index: 2;
    transform-origin: center bottom;
    display: flex;
    flex-direction: column;
    align-items: center;
  }

  .flask__neck {
    width: 1.2em;
    height: 1.3em;
    margin-bottom: -0.3em;
    position: relative;
    z-index: 1;
    background: linear-gradient(
      to right,
      rgba(190, 220, 240, 0.5),
      rgba(255, 255, 255, 0.85) 25% 75%,
      rgba(190, 220, 240, 0.5)
    );
    border: 0.09em solid var(--glass-edge);
    border-bottom: none;
    border-radius: 0.18em 0.18em 0 0;
    box-shadow:
      inset 0.08em 0.08em 0.15em rgba(255, 255, 255, 0.7),
      inset -0.05em -0.05em 0.12em rgba(0, 0, 0, 0.04),
      0 -0.05em 0.15em rgba(0, 0, 0, 0.06);
  }

  .flask__neck::before {
    content: "";
    position: absolute;
    top: -0.2em;
    left: 50%;
    transform: translateX(-50%);
    width: 1.65em;
    height: 0.28em;
    background: linear-gradient(
      to bottom,
      rgba(240, 248, 255, 0.95),
      rgba(200, 225, 245, 0.8)
    );
    border-radius: 0.14em;
    border: 0.06em solid var(--glass-edge);
    box-shadow:
      inset 0.04em 0.04em 0.08em rgba(255, 255, 255, 0.9),
      inset -0.03em -0.03em 0.06em rgba(0, 0, 0, 0.06),
      0 0.06em 0.12em rgba(0, 0, 0, 0.08);
  }

  .flask__neck::after {
    content: "";
    position: absolute;
    bottom: 0;
    left: 0.12em;
    right: 0.12em;
    height: 28%;
    background: var(--liquid-day-surface);
    border-radius: 0 0 0.04em 0.04em;
    transition: var(--transition);
  }

  /* ── Vapor ── */

  .flask__vapor {
    position: absolute;
    border-radius: var(--radius);
    transition: var(--transition);
    opacity: 0;
  }

  .flask__vapor:nth-child(1) {
    width: 0.4em;
    height: 0.4em;
    top: -0.5em;
    left: 0.1em;
    background: radial-gradient(circle, rgba(255, 255, 255, 0.45), rgba(255, 255, 255, 0.1) 70%, transparent);
  }

  .flask__vapor:nth-child(2) {
    width: 0.3em;
    height: 0.3em;
    top: -0.85em;
    right: 0.15em;
    background: radial-gradient(circle, rgba(255, 255, 255, 0.35), rgba(255, 255, 255, 0.08) 70%, transparent);
  }

  /* ── Flask body (claymorphism glass sphere) ── */

  .flask__body {
    width: var(--flask-diameter);
    height: var(--flask-diameter);
    border-radius: var(--radius);
    position: relative;
    overflow: hidden;
    transition: var(--transition);
    z-index: 1;
    background: radial-gradient(
      ellipse at 32% 28%,
      rgba(255, 255, 255, 0.98) 0%,
      rgba(240, 248, 255, 0.9) 20%,
      rgba(220, 238, 252, 0.82) 45%,
      rgba(195, 220, 240, 0.75) 70%,
      rgba(175, 205, 230, 0.7) 100%
    );
    border: 0.1em solid rgba(180, 210, 235, 0.35);
    box-shadow:
      inset 0.25em 0.25em 0.6em rgba(255, 255, 255, 0.7),
      inset -0.2em -0.2em 0.5em rgba(100, 150, 200, 0.12),
      inset 0 0.15em 0.3em rgba(255, 255, 255, 0.4),
      0 0.25em 0.6em rgba(0, 0, 0, 0.1),
      0 0.08em 0.2em rgba(0, 0, 0, 0.06);
  }

  .flask__body::before {
    content: "";
    position: absolute;
    bottom: 0;
    left: 0;
    width: 100%;
    height: 56%;
    background:
      radial-gradient(
        ellipse at 50% 0%,
        rgba(255, 255, 255, 0.3) 0%,
        transparent 60%
      ),
      linear-gradient(
        to bottom,
        var(--liquid-day-surface) 0%,
        rgba(13, 148, 136, 0.65) 20%,
        var(--liquid-day) 60%,
        rgba(6, 95, 86, 0.95) 100%
      );
    border-radius: 0 0 var(--radius) var(--radius);
    transition: var(--transition);
    box-shadow:
      inset 0 0.2em 0.4em rgba(255, 255, 255, 0.15),
      inset 0 -0.3em 0.6em rgba(0, 60, 50, 0.2);
  }

  .flask__body::after {
    content: "";
    position: absolute;
    top: 0.4em;
    left: 0.5em;
    width: 1.1em;
    height: 0.65em;
    background: radial-gradient(
      ellipse,
      rgba(255, 255, 255, 0.7) 0%,
      rgba(255, 255, 255, 0.3) 50%,
      transparent 100%
    );
    border-radius: 50%;
    transform: rotate(-20deg);
    transition: var(--transition);
    box-shadow:
      1.1em 2.3em 0 -0.04em rgba(255, 255, 255, 0.5),
      2.2em 2.5em 0 -0.05em rgba(255, 255, 255, 0.4),
      0.5em 2.7em 0 -0.03em rgba(255, 255, 255, 0.35),
      1.7em 2.0em 0 -0.07em rgba(255, 255, 255, 0.3),
      2.6em 2.2em 0 -0.06em rgba(255, 255, 255, 0.35),
      0.8em 1.9em 0 -0.08em rgba(255, 255, 255, 0.25),
      1.4em 2.9em 0 -0.09em rgba(255, 255, 255, 0.3),
      2.0em 3.0em 0 -0.1em rgba(255, 255, 255, 0.2),
      2.2em 0.6em 0 -0.12em rgba(255, 255, 255, 0.2);
  }

  .artificial__hidden {
    position: absolute;
    border-radius: inherit;
    inset: 0;
    pointer-events: none;
    overflow: hidden;
  }

  .flask__shadow {
    width: var(--flask-diameter);
    height: 20%;
    border-radius: 50%;
    background: #2a3a4c;
    box-shadow: 0.313em 0 3.125em #2a3a4c;
    opacity: 0.25;
    position: absolute;
    bottom: 0;
    left: calc(var(--toggle-offset) - 0.938em);
    transition: var(--transition);
    transform: skew(-70deg);
    z-index: 1;
  }

  /* ── Scenery ── */

  .chem-toggle__scenery {
    width: 100%;
    height: 100%;
    pointer-events: none;
    overflow: hidden;
    position: relative;
    border-radius: inherit;
  }

  /* Ground */
  .chem-toggle__scenery::before {
    content: "";
    position: absolute;
    width: 100%;
    height: 30%;
    bottom: 0;
    background: linear-gradient(to bottom, #8ba4b8, #6b8a9e);
    z-index: 1;
  }

  /* Clouds */
  .chem-toggle__cloud {
    z-index: 1;
    position: absolute;
    border-radius: 50%;
  }

  .chem-toggle__cloud:nth-last-child(1) {
    width: 0.875em;
    height: 0.625em;
    filter: blur(0.125em) drop-shadow(0.313em 0.313em #ffffffae)
      drop-shadow(-0.625em 0 #fff) drop-shadow(-0.938em -0.125em #fff);
    right: 1.875em;
    top: 2.813em;
    background: linear-gradient(to top right, #ffffffae, #ffffffae);
    transition: var(--transition);
  }

  .chem-toggle__cloud:nth-last-child(2) {
    top: 0.625em;
    right: 4.375em;
    width: 0.875em;
    height: 0.375em;
    background: #dfdedeae;
    filter: blur(0.125em) drop-shadow(-0.313em -0.188em #e0dfdfae)
      drop-shadow(-0.625em -0.188em #bbbbbbae) drop-shadow(-1em 0.063em #cfcfcfae);
    transition: 0.6s;
  }

  .chem-toggle__cloud:nth-last-child(3) {
    top: 1.25em;
    right: 0.938em;
    width: 0.875em;
    height: 0.375em;
    background: #ffffffae;
    filter: blur(0.125em) drop-shadow(0.438em 0.188em #ffffffae)
      drop-shadow(-0.625em 0.313em #ffffffae);
    transition: 0.8s;
  }

  /* Moon + craters */
  .moon,
  .moon-crater-1,
  .moon-crater-2 {
    position: absolute;
    border-radius: var(--radius);
    background: linear-gradient(#fff, #6e8ea2);
    top: 100%;
  }

  .moon {
    left: 0.938em;
    width: 1.875em;
    height: 1.875em;
    box-shadow: 0 0 0.188em #ffffff52, 0 0 0.188em #6e8ea24b;
    transition: var(--transition);
  }

  .moon::before,
  .moon::after {
    content: "";
    position: absolute;
    border-radius: inherit;
    box-shadow: inset 0 0 0.063em rgb(140, 162, 169);
    background: rgb(184, 196, 200);
  }

  .moon::before {
    left: 0.313em;
    top: 0.313em;
    width: 0.438em;
    height: 0.438em;
  }

  .moon::after {
    width: 0.25em;
    height: 0.25em;
    left: 1.25em;
    top: 0.75em;
  }

  .moon-crater-1,
  .moon-crater-2 {
    box-shadow: 0 0 0.125em #ffffff52, 0 0 0.125em #6e8ea24b;
  }

  .moon-crater-1 {
    left: 3.438em;
    width: 0.625em;
    height: 0.625em;
    transition: 0.6s;
  }

  .moon-crater-2 {
    left: 4.375em;
    width: 0.5em;
    height: 0.5em;
    transition: 0.8s;
  }

  /* Sun */
  .sun-primary,
  .sun-secondary {
    position: absolute;
    width: 1.25em;
    height: 1.25em;
    border-radius: var(--radius);
  }

  .sun-primary {
    background: #fefefe;
    right: 3.125em;
    top: 0.625em;
    box-shadow: 0 0 0.438em #fdf4e1;
    transition: var(--transition);
  }

  .sun-secondary {
    background: linear-gradient(#e6ac5c, #d75449);
    right: 1.25em;
    top: 2.188em;
    box-shadow: 0 0 0.438em #e6ad5c3d, 0 0 0.438em #d755494f;
    transition: 0.7s;
  }

  /* Stars */
  .chem-toggle__star {
    position: absolute;
    width: 0.063em;
    height: 0.063em;
    background: #fff;
    border-radius: var(--radius);
    filter: drop-shadow(0 0 0.063em #fff);
    color: #fff;
    top: 100%;
  }

  .chem-toggle__star:nth-child(1) {
    left: 3.75em;
    box-shadow: 1.25em 0.938em, -1.25em 2.5em, 0 1.25em, 1.875em 0.625em,
      -3.125em 1.875em, 1.25em 2.813em;
    transition: 0.2s;
  }

  .chem-toggle__star:nth-child(2) {
    left: 4.688em;
    box-shadow: 0.625em 0, 0 0.625em, -0.625em -0.625em, 0.625em 0.938em,
      -3.125em 1.25em, 1.25em -1.563em;
    transition: 0.3s;
  }

  .chem-toggle__star:nth-child(3) {
    left: 5.313em;
    box-shadow: -0.625em -0.625em, -2.188em 1.25em, -2.188em 0, -3.75em -0.625em,
      -3.125em -0.625em, -2.5em -0.313em, 0.75em -0.625em;
    transition: var(--transition);
  }

  .chem-toggle__star:nth-child(4),
  .chem-toggle__star:nth-child(5),
  .chem-toggle__star:nth-child(6),
  .chem-toggle__star:nth-child(7) {
    width: 0.125em;
    height: 0.125em;
  }

  .chem-toggle__star:nth-child(4) {
    left: 1.875em;
    transition: 0.5s;
  }

  .chem-toggle__star:nth-child(5) {
    left: 5em;
    transition: 0.6s;
  }

  .chem-toggle__star:nth-child(6) {
    left: 2.5em;
    transition: 0.7s;
  }

  .chem-toggle__star:nth-child(7) {
    left: 3.438em;
    transition: 0.8s;
  }

  /* ── Checked (dark mode) state ── */

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .chem-toggle__star:nth-child(1) {
    top: 0.625em;
  }

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .chem-toggle__star:nth-child(2) {
    top: 1.875em;
  }

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .chem-toggle__star:nth-child(3) {
    top: 1.25em;
  }

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .chem-toggle__star:nth-child(4) {
    top: 3.438em;
  }

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .chem-toggle__star:nth-child(5) {
    top: 3.438em;
  }

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .chem-toggle__star:nth-child(6) {
    top: 0.313em;
  }

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .chem-toggle__star:nth-child(7) {
    top: 1.875em;
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .chem-toggle__cloud {
    right: -100%;
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .moon {
    top: 0.938em;
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .moon-crater-1 {
    top: 2.5em;
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .moon-crater-2 {
    top: 2.75em;
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container {
    background-position-y: 0;
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .sun-primary,
  .chem-toggle__checkbox:checked + .chem-toggle__container .sun-secondary {
    top: 100%;
  }

  /* Flask slides right */
  .chem-toggle__checkbox:checked + .chem-toggle__container .flask {
    left: calc(100% - var(--flask-diameter) - var(--toggle-offset));
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .flask__shadow {
    left: calc(100% - var(--flask-diameter) - var(--toggle-offset) + 0.938em);
    transform: skew(70deg);
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .flask__body {
    transform: rotate(360deg);
    background: radial-gradient(
      ellipse at 32% 28%,
      rgba(230, 240, 255, 0.92) 0%,
      rgba(210, 225, 245, 0.82) 20%,
      rgba(185, 210, 235, 0.72) 45%,
      rgba(160, 190, 220, 0.65) 70%,
      rgba(140, 175, 210, 0.6) 100%
    );
    box-shadow:
      inset 0.25em 0.25em 0.6em rgba(200, 220, 255, 0.5),
      inset -0.2em -0.2em 0.5em rgba(50, 80, 120, 0.15),
      inset 0 0.15em 0.3em rgba(200, 220, 255, 0.3),
      0 0.25em 0.6em rgba(0, 0, 0, 0.15),
      0 0.08em 0.2em rgba(0, 0, 0, 0.1);
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .flask__body::before {
    background:
      radial-gradient(
        ellipse at 50% 0%,
        rgba(100, 140, 180, 0.2) 0%,
        transparent 55%
      ),
      linear-gradient(
        to bottom,
        var(--liquid-night-surface) 0%,
        rgba(30, 41, 59, 0.75) 20%,
        var(--liquid-night) 55%,
        var(--liquid-night-deep) 100%
      );
    box-shadow:
      inset 0 0.15em 0.3em rgba(100, 140, 200, 0.1),
      inset 0 -0.3em 0.6em rgba(0, 0, 0, 0.25);
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .flask__body::after {
    box-shadow:
      1.1em 2.3em 0 -0.04em rgba(150, 180, 220, 0.3),
      2.2em 2.5em 0 -0.05em rgba(150, 180, 220, 0.25),
      0.5em 2.7em 0 -0.03em rgba(150, 180, 220, 0.2),
      1.7em 2.0em 0 -0.07em rgba(150, 180, 220, 0.18),
      2.6em 2.2em 0 -0.06em rgba(150, 180, 220, 0.2),
      0.8em 1.9em 0 -0.08em rgba(150, 180, 220, 0.15),
      1.4em 2.9em 0 -0.09em rgba(150, 180, 220, 0.18),
      2.0em 3.0em 0 -0.1em rgba(150, 180, 220, 0.12),
      2.2em 0.6em 0 -0.12em rgba(180, 200, 230, 0.15);
  }

  .chem-toggle__checkbox:checked
    + .chem-toggle__container
    .flask__neck::after {
    background: var(--liquid-night-surface);
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container .flask__vapor {
    opacity: 0.35;
  }

  /* ── Hover interactions ── */

  .chem-toggle:hover .flask__body::after {
    box-shadow:
      1.1em 2.0em 0 -0.03em rgba(255, 255, 255, 0.55),
      2.2em 2.2em 0 -0.04em rgba(255, 255, 255, 0.45),
      0.5em 2.4em 0 -0.02em rgba(255, 255, 255, 0.4),
      1.7em 1.7em 0 -0.06em rgba(255, 255, 255, 0.35),
      2.6em 1.9em 0 -0.05em rgba(255, 255, 255, 0.4),
      0.8em 1.6em 0 -0.07em rgba(255, 255, 255, 0.3),
      1.4em 2.6em 0 -0.08em rgba(255, 255, 255, 0.35),
      2.0em 2.7em 0 -0.09em rgba(255, 255, 255, 0.25),
      2.2em 0.6em 0 -0.12em rgba(255, 255, 255, 0.25);
  }

  .chem-toggle__checkbox:checked + .chem-toggle__container:hover .flask__body::after {
    box-shadow:
      1.1em 2.0em 0 -0.03em rgba(150, 180, 220, 0.35),
      2.2em 2.2em 0 -0.04em rgba(150, 180, 220, 0.3),
      0.5em 2.4em 0 -0.02em rgba(150, 180, 220, 0.25),
      1.7em 1.7em 0 -0.06em rgba(150, 180, 220, 0.22),
      2.6em 1.9em 0 -0.05em rgba(150, 180, 220, 0.25),
      0.8em 1.6em 0 -0.07em rgba(150, 180, 220, 0.18),
      1.4em 2.6em 0 -0.08em rgba(150, 180, 220, 0.22),
      2.0em 2.7em 0 -0.09em rgba(150, 180, 220, 0.15),
      2.2em 0.6em 0 -0.12em rgba(180, 200, 230, 0.18);
  }

  .chem-toggle:hover .flask__vapor {
    opacity: 0.45;
  }

  .chem-toggle:hover .flask__vapor:nth-child(1) {
    top: -0.75em;
    opacity: 0.4;
  }

  .chem-toggle:hover .flask__vapor:nth-child(2) {
    top: -1.15em;
    opacity: 0.3;
  }

  .chem-toggle__checkbox:active
    + .chem-toggle__container
    .flask__neck-container {
    transform: rotate(18deg);
  }

  .chem-toggle__checkbox:checked:active
    + .chem-toggle__container
    .flask__neck-container {
    transform: rotate(-18deg);
  }

  .chem-toggle__checkbox:active + .chem-toggle__container .flask__body {
    box-shadow:
      inset 0.3em 0.3em 0.7em rgba(255, 255, 255, 0.5),
      inset -0.25em -0.25em 0.6em rgba(100, 150, 200, 0.15),
      inset 0 0.2em 0.35em rgba(255, 255, 255, 0.3),
      0 0.15em 0.35em rgba(0, 0, 0, 0.12),
      0 0.04em 0.1em rgba(0, 0, 0, 0.08);
  }
`;
