/**
 * BB8ThemeToggle
 *
 * Star Wars BB-8 droid toggle switch for light/dark mode.
 * Unchecked = light (day scene), checked = dark (night scene).
 * Wrapped in a claymorphism container to match the app design system.
 */

import styled from 'styled-components';
import { useThemeContext } from '../../contexts/ThemeContext';
import { cn } from '../../lib/utils';

interface BB8ThemeToggleProps {
  className?: string;
}

export function BB8ThemeToggle({ className }: BB8ThemeToggleProps) {
  const { resolvedTheme, setTheme } = useThemeContext();
  const isDark = resolvedTheme === 'dark';

  const handleChange = () => {
    setTheme(isDark ? 'light' : 'dark');
  };

  return (
    <StyledWrapper
      className={cn(
        'overflow-visible',
        'transition-all duration-300',
        className,
      )}
    >
      <label className="bb8-toggle" aria-label={`Switch to ${isDark ? 'light' : 'dark'} mode`}>
        <input
          className="bb8-toggle__checkbox"
          type="checkbox"
          checked={isDark}
          onChange={handleChange}
        />
        <div className="bb8-toggle__container">
          <div className="bb8-toggle__scenery">
            <div className="bb8-toggle__star" />
            <div className="bb8-toggle__star" />
            <div className="bb8-toggle__star" />
            <div className="bb8-toggle__star" />
            <div className="bb8-toggle__star" />
            <div className="bb8-toggle__star" />
            <div className="bb8-toggle__star" />
            <div className="tatto-1" />
            <div className="tatto-2" />
            <div className="gomrassen" />
            <div className="hermes" />
            <div className="chenini" />
            <div className="bb8-toggle__cloud" />
            <div className="bb8-toggle__cloud" />
            <div className="bb8-toggle__cloud" />
          </div>
          <div className="bb8">
            <div className="bb8__head-container">
              <div className="bb8__antenna" />
              <div className="bb8__antenna" />
              <div className="bb8__head" />
            </div>
            <div className="bb8__body" />
          </div>
          <div className="artificial__hidden">
            <div className="bb8__shadow" />
          </div>
        </div>
      </label>
    </StyledWrapper>
  );
}

const StyledWrapper = styled.div`
  .bb8-toggle {
    --toggle-size: 5.25px;
    --toggle-width: 10.625em;
    --toggle-height: 5.625em;
    --toggle-offset: calc((var(--toggle-height) - var(--bb8-diameter)) / 2);
    --toggle-bg: linear-gradient(#2c4770, #070e2b 35%, #628cac 50% 70%, #a6c5d4)
      no-repeat;
    --bb8-diameter: 4.375em;
    --radius: 99em;
    --transition: 0.4s;
    --accent: #de7d2f;
    --bb8-bg: #fff;
  }

  .bb8-toggle,
  .bb8-toggle *,
  .bb8-toggle *::before,
  .bb8-toggle *::after {
    box-sizing: border-box;
  }

  .bb8-toggle {
    cursor: pointer;
    font-size: var(--toggle-size);
  }

  .bb8-toggle__checkbox {
    display: none;
  }

  .bb8-toggle__container {
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

  .bb8 {
    display: flex;
    flex-direction: column;
    align-items: center;
    position: absolute;
    top: calc(var(--toggle-offset) - 1.688em + 0.188em);
    left: var(--toggle-offset);
    transition: var(--transition);
    z-index: 2;
  }

  .bb8__head-container {
    position: relative;
    transition: var(--transition);
    z-index: 2;
    transform-origin: 1.25em 3.75em;
  }

  .bb8__head {
    overflow: hidden;
    margin-bottom: -0.188em;
    width: 2.5em;
    height: 1.688em;
    background: linear-gradient(
        transparent 0.063em,
        dimgray 0.063em 0.313em,
        transparent 0.313em 0.375em,
        var(--accent) 0.375em 0.5em,
        transparent 0.5em 1.313em,
        silver 1.313em 1.438em,
        transparent 1.438em
      ),
      linear-gradient(
        45deg,
        transparent 0.188em,
        var(--bb8-bg) 0.188em 1.25em,
        transparent 1.25em
      ),
      linear-gradient(
        -45deg,
        transparent 0.188em,
        var(--bb8-bg) 0.188em 1.25em,
        transparent 1.25em
      ),
      linear-gradient(var(--bb8-bg) 1.25em, transparent 1.25em);
    border-radius: var(--radius) var(--radius) 0 0;
    position: relative;
    z-index: 1;
    filter: drop-shadow(0 0.063em 0.125em gray);
  }

  .bb8__head::before {
    content: "";
    position: absolute;
    width: 0.563em;
    height: 0.563em;
    background: radial-gradient(
        0.125em circle at 0.25em 0.375em,
        red,
        transparent
      ),
      radial-gradient(
        0.063em circle at 0.375em 0.188em,
        var(--bb8-bg) 50%,
        transparent 100%
      ),
      linear-gradient(45deg, #000 0.188em, dimgray 0.313em 0.375em, #000 0.5em);
    border-radius: var(--radius);
    top: 0.413em;
    left: 50%;
    transform: translate(-50%);
    box-shadow: 0 0 0 0.089em lightgray, 0.563em 0.281em 0 -0.148em,
      0.563em 0.281em 0 -0.1em var(--bb8-bg), 0.563em 0.281em 0 -0.063em;
    z-index: 1;
    transition: var(--transition);
  }

  .bb8__head::after {
    content: "";
    position: absolute;
    bottom: 0.375em;
    left: 0;
    width: 100%;
    height: 0.188em;
    background: linear-gradient(
      to right,
      var(--accent) 0.125em,
      transparent 0.125em 0.188em,
      var(--accent) 0.188em 0.313em,
      transparent 0.313em 0.375em,
      var(--accent) 0.375em 0.938em,
      transparent 0.938em 1em,
      var(--accent) 1em 1.125em,
      transparent 1.125em 1.875em,
      var(--accent) 1.875em 2em,
      transparent 2em 2.063em,
      var(--accent) 2.063em 2.25em,
      transparent 2.25em 2.313em,
      var(--accent) 2.313em 2.375em,
      transparent 2.375em 2.438em,
      var(--accent) 2.438em
    );
    transition: var(--transition);
  }

  .bb8__antenna {
    position: absolute;
    transform: translateY(-90%);
    width: 0.059em;
    border-radius: var(--radius) var(--radius) 0 0;
    transition: var(--transition);
  }

  .bb8__antenna:nth-child(1) {
    height: 0.938em;
    right: 0.938em;
    background: linear-gradient(#000 0.188em, silver 0.188em);
  }

  .bb8__antenna:nth-child(2) {
    height: 0.375em;
    left: 50%;
    transform: translate(-50%, -90%);
    background: silver;
  }

  .bb8__body {
    width: 4.375em;
    height: 4.375em;
    border-radius: var(--radius);
    position: relative;
    overflow: hidden;
    transition: var(--transition);
    z-index: 1;
    transform: rotate(45deg);
    background: linear-gradient(
        -90deg,
        var(--bb8-bg) 4%,
        var(--accent) 4% 10%,
        transparent 10% 90%,
        var(--accent) 90% 96%,
        var(--bb8-bg) 96%
      ),
      linear-gradient(
        var(--bb8-bg) 4%,
        var(--accent) 4% 10%,
        transparent 10% 90%,
        var(--accent) 90% 96%,
        var(--bb8-bg) 96%
      ),
      linear-gradient(
        to right,
        transparent 2.156em,
        silver 2.156em 2.219em,
        transparent 2.188em
      ),
      linear-gradient(
        transparent 2.156em,
        silver 2.156em 2.219em,
        transparent 2.188em
      );
    background-color: var(--bb8-bg);
  }

  .bb8__body::after {
    content: "";
    bottom: 1.5em;
    left: 0.563em;
    position: absolute;
    width: 0.188em;
    height: 0.188em;
    background: rgb(236, 236, 236);
    color: rgb(236, 236, 236);
    border-radius: 50%;
    box-shadow: 0.875em 0.938em, 0 -1.25em, 0.875em -2.125em, 2.125em -2.125em,
      3.063em -1.25em, 3.063em 0, 2.125em 0.938em;
  }

  .bb8__body::before {
    content: "";
    width: 2.625em;
    height: 2.625em;
    position: absolute;
    border-radius: 50%;
    z-index: 0.1;
    overflow: hidden;
    top: 50%;
    left: 50%;
    transform: translate(-50%, -50%);
    border: 0.313em solid var(--accent);
    background: radial-gradient(
        1em circle at center,
        rgb(236, 236, 236) 50%,
        transparent 51%
      ),
      radial-gradient(1.25em circle at center, var(--bb8-bg) 50%, transparent 51%),
      linear-gradient(
        -90deg,
        transparent 42%,
        var(--accent) 42% 58%,
        transparent 58%
      ),
      linear-gradient(var(--bb8-bg) 42%, var(--accent) 42% 58%, var(--bb8-bg) 58%);
  }

  .artificial__hidden {
    position: absolute;
    border-radius: inherit;
    inset: 0;
    pointer-events: none;
    overflow: hidden;
  }

  .bb8__shadow {
    width: var(--bb8-diameter);
    height: 20%;
    border-radius: 50%;
    background: #3a271c;
    box-shadow: 0.313em 0 3.125em #3a271c;
    opacity: 0.25;
    position: absolute;
    bottom: 0;
    left: calc(var(--toggle-offset) - 0.938em);
    transition: var(--transition);
    transform: skew(-70deg);
    z-index: 1;
  }

  .bb8-toggle__scenery {
    width: 100%;
    height: 100%;
    pointer-events: none;
    overflow: hidden;
    position: relative;
    border-radius: inherit;
  }

  .bb8-toggle__scenery::before {
    content: "";
    position: absolute;
    width: 100%;
    height: 30%;
    bottom: 0;
    background: #b18d71;
    z-index: 1;
  }

  .bb8-toggle__cloud {
    z-index: 1;
    position: absolute;
    border-radius: 50%;
  }

  .bb8-toggle__cloud:nth-last-child(1) {
    width: 0.875em;
    height: 0.625em;
    filter: blur(0.125em) drop-shadow(0.313em 0.313em #ffffffae)
      drop-shadow(-0.625em 0 #fff) drop-shadow(-0.938em -0.125em #fff);
    right: 1.875em;
    top: 2.813em;
    background: linear-gradient(to top right, #ffffffae, #ffffffae);
    transition: var(--transition);
  }

  .bb8-toggle__cloud:nth-last-child(2) {
    top: 0.625em;
    right: 4.375em;
    width: 0.875em;
    height: 0.375em;
    background: #dfdedeae;
    filter: blur(0.125em) drop-shadow(-0.313em -0.188em #e0dfdfae)
      drop-shadow(-0.625em -0.188em #bbbbbbae) drop-shadow(-1em 0.063em #cfcfcfae);
    transition: 0.6s;
  }

  .bb8-toggle__cloud:nth-last-child(3) {
    top: 1.25em;
    right: 0.938em;
    width: 0.875em;
    height: 0.375em;
    background: #ffffffae;
    filter: blur(0.125em) drop-shadow(0.438em 0.188em #ffffffae)
      drop-shadow(-0.625em 0.313em #ffffffae);
    transition: 0.8s;
  }

  .gomrassen,
  .hermes,
  .chenini {
    position: absolute;
    border-radius: var(--radius);
    background: linear-gradient(#fff, #6e8ea2);
    top: 100%;
  }

  .gomrassen {
    left: 0.938em;
    width: 1.875em;
    height: 1.875em;
    box-shadow: 0 0 0.188em #ffffff52, 0 0 0.188em #6e8ea24b;
    transition: var(--transition);
  }

  .gomrassen::before,
  .gomrassen::after {
    content: "";
    position: absolute;
    border-radius: inherit;
    box-shadow: inset 0 0 0.063em rgb(140, 162, 169);
    background: rgb(184, 196, 200);
  }

  .gomrassen::before {
    left: 0.313em;
    top: 0.313em;
    width: 0.438em;
    height: 0.438em;
  }

  .gomrassen::after {
    width: 0.25em;
    height: 0.25em;
    left: 1.25em;
    top: 0.75em;
  }

  .hermes {
    left: 3.438em;
    width: 0.625em;
    height: 0.625em;
    box-shadow: 0 0 0.125em #ffffff52, 0 0 0.125em #6e8ea24b;
    transition: 0.6s;
  }

  .chenini {
    left: 4.375em;
    width: 0.5em;
    height: 0.5em;
    box-shadow: 0 0 0.125em #ffffff52, 0 0 0.125em #6e8ea24b;
    transition: 0.8s;
  }

  .tatto-1,
  .tatto-2 {
    position: absolute;
    width: 1.25em;
    height: 1.25em;
    border-radius: var(--radius);
  }

  .tatto-1 {
    background: #fefefe;
    right: 3.125em;
    top: 0.625em;
    box-shadow: 0 0 0.438em #fdf4e1;
    transition: var(--transition);
  }

  .tatto-2 {
    background: linear-gradient(#e6ac5c, #d75449);
    right: 1.25em;
    top: 2.188em;
    box-shadow: 0 0 0.438em #e6ad5c3d, 0 0 0.438em #d755494f;
    transition: 0.7s;
  }

  .bb8-toggle__star {
    position: absolute;
    width: 0.063em;
    height: 0.063em;
    background: #fff;
    border-radius: var(--radius);
    filter: drop-shadow(0 0 0.063em #fff);
    color: #fff;
    top: 100%;
  }

  .bb8-toggle__star:nth-child(1) {
    left: 3.75em;
    box-shadow: 1.25em 0.938em, -1.25em 2.5em, 0 1.25em, 1.875em 0.625em,
      -3.125em 1.875em, 1.25em 2.813em;
    transition: 0.2s;
  }

  .bb8-toggle__star:nth-child(2) {
    left: 4.688em;
    box-shadow: 0.625em 0, 0 0.625em, -0.625em -0.625em, 0.625em 0.938em,
      -3.125em 1.25em, 1.25em -1.563em;
    transition: 0.3s;
  }

  .bb8-toggle__star:nth-child(3) {
    left: 5.313em;
    box-shadow: -0.625em -0.625em, -2.188em 1.25em, -2.188em 0, -3.75em -0.625em,
      -3.125em -0.625em, -2.5em -0.313em, 0.75em -0.625em;
    transition: var(--transition);
  }

  .bb8-toggle__star:nth-child(4) {
    left: 1.875em;
    width: 0.125em;
    height: 0.125em;
    transition: 0.5s;
  }

  .bb8-toggle__star:nth-child(5) {
    left: 5em;
    width: 0.125em;
    height: 0.125em;
    transition: 0.6s;
  }

  .bb8-toggle__star:nth-child(6) {
    left: 2.5em;
    width: 0.125em;
    height: 0.125em;
    transition: 0.7s;
  }

  .bb8-toggle__star:nth-child(7) {
    left: 3.438em;
    width: 0.125em;
    height: 0.125em;
    transition: 0.8s;
  }

  /* ── Checked (dark mode) state ── */

  .bb8-toggle__checkbox:checked
    + .bb8-toggle__container
    .bb8-toggle__star:nth-child(1) {
    top: 0.625em;
  }

  .bb8-toggle__checkbox:checked
    + .bb8-toggle__container
    .bb8-toggle__star:nth-child(2) {
    top: 1.875em;
  }

  .bb8-toggle__checkbox:checked
    + .bb8-toggle__container
    .bb8-toggle__star:nth-child(3) {
    top: 1.25em;
  }

  .bb8-toggle__checkbox:checked
    + .bb8-toggle__container
    .bb8-toggle__star:nth-child(4) {
    top: 3.438em;
  }

  .bb8-toggle__checkbox:checked
    + .bb8-toggle__container
    .bb8-toggle__star:nth-child(5) {
    top: 3.438em;
  }

  .bb8-toggle__checkbox:checked
    + .bb8-toggle__container
    .bb8-toggle__star:nth-child(6) {
    top: 0.313em;
  }

  .bb8-toggle__checkbox:checked
    + .bb8-toggle__container
    .bb8-toggle__star:nth-child(7) {
    top: 1.875em;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .bb8-toggle__cloud {
    right: -100%;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .gomrassen {
    top: 0.938em;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .hermes {
    top: 2.5em;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .chenini {
    top: 2.75em;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container {
    background-position-y: 0;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .tatto-1 {
    top: 100%;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .tatto-2 {
    top: 100%;
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .bb8 {
    left: calc(100% - var(--bb8-diameter) - var(--toggle-offset));
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .bb8__shadow {
    left: calc(100% - var(--bb8-diameter) - var(--toggle-offset) + 0.938em);
    transform: skew(70deg);
  }

  .bb8-toggle__checkbox:checked + .bb8-toggle__container .bb8__body {
    transform: rotate(225deg);
  }

  /* ── Hover interactions ── */

  .bb8-toggle__checkbox:hover + .bb8-toggle__container .bb8__head::before {
    left: 100%;
  }

  .bb8-toggle__checkbox:not(:checked):hover
    + .bb8-toggle__container
    .bb8__antenna:nth-child(1) {
    right: 1.5em;
  }

  .bb8-toggle__checkbox:hover
    + .bb8-toggle__container
    .bb8__antenna:nth-child(2) {
    left: 0.938em;
  }

  .bb8-toggle__checkbox:hover + .bb8-toggle__container .bb8__head::after {
    background-position: 1.375em 0;
  }

  .bb8-toggle__checkbox:checked:hover
    + .bb8-toggle__container
    .bb8__head::before {
    left: 0;
  }

  .bb8-toggle__checkbox:checked:hover
    + .bb8-toggle__container
    .bb8__antenna:nth-child(2) {
    left: calc(100% - 0.938em);
  }

  .bb8-toggle__checkbox:checked:hover
    + .bb8-toggle__container
    .bb8__head::after {
    background-position: -1.375em 0;
  }

  .bb8-toggle__checkbox:active + .bb8-toggle__container .bb8__head-container {
    transform: rotate(25deg);
  }

  .bb8-toggle__checkbox:checked:active
    + .bb8-toggle__container
    .bb8__head-container {
    transform: rotate(-25deg);
  }

  .bb8:hover .bb8__head::before,
  .bb8:hover .bb8__antenna:nth-child(2) {
    left: 50% !important;
  }

  .bb8:hover .bb8__antenna:nth-child(1) {
    right: 0.938em !important;
  }

  .bb8:hover .bb8__head::after {
    background-position: 0 0 !important;
  }
`;
