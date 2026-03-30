import { useState, useEffect, useRef, useCallback } from "react";

// ── PARAMETERS ──────────────────────────────────────────────────────────────
const CELL_RADIUS       = 10;      // px
const SPRING_K          = 18;
const DAMPING           = 1.0;
const CUTOFF            = 1.5 * 2 * CELL_RADIUS;
const DT                = 0.08;
const MEAN_CYCLE        = 280;     // frames
const SD_CYCLE          = 30;
const BIRTH_SEP         = 2;
const MAX_CELLS         = 120;

// Generation colour palette — vivid biological hues
const GEN_COLORS = [
  "#00e5ff", "#69ff47", "#ff6b6b", "#ffd93d",
  "#c77dff", "#ff9a3c", "#06d6a0", "#f72585",
  "#4cc9f0", "#b5e48c"
];

let _idCounter = 0;

function randn() {
  // Box-Muller
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function makeCycle(birthFrame) {
  return {
    birthFrame,
    duration: Math.max(120, MEAN_CYCLE + randn() * SD_CYCLE),
  };
}

function makeCell(x, y, generation, birthFrame) {
  return {
    id: _idCounter++,
    x, y,
    vx: 0, vy: 0,
    generation,
    cycle: makeCycle(birthFrame),
    dividing: false,    // flash state
    flashTimer: 0,
    // For smooth birth scale-in
    scale: 0.3,
  };
}

function computeForces(cells) {
  const fx = new Float64Array(cells.length);
  const fy = new Float64Array(cells.length);
  for (let i = 0; i < cells.length; i++) {
    for (let j = i + 1; j < cells.length; j++) {
      const dx = cells[j].x - cells[i].x;
      const dy = cells[j].y - cells[i].y;
      const dist = Math.sqrt(dx * dx + dy * dy) || 0.001;
      if (dist > CUTOFF) continue;
      const nat = CELL_RADIUS * 2;
      const overlap = dist - nat;
      const f = SPRING_K * overlap / dist;
      fx[i] += f * dx; fy[i] += f * dy;
      fx[j] -= f * dx; fy[j] -= f * dy;
    }
  }
  return { fx, fy };
}

export default function ChaseSim() {
  const canvasRef   = useRef(null);
  const stateRef    = useRef({
    cells: [],
    frame: 0,
    running: false,
    speed: 1,
    width: 600,
    height: 600,
  });
  const rafRef      = useRef(null);
  const [uiState, setUiState] = useState({
    frame: 0,
    count: 1,
    running: false,
    speed: 1,
    maxReached: false,
  });

  // ── Init ──────────────────────────────────────────────────────────────────
  const reset = useCallback(() => {
    _idCounter = 0;
    const s = stateRef.current;
    s.cells = [makeCell(s.width / 2, s.height / 2, 0, 0)];
    s.frame = 0;
    setUiState(u => ({ ...u, frame: 0, count: 1, maxReached: false }));
  }, []);

  // ── Simulation step ───────────────────────────────────────────────────────
  const step = useCallback(() => {
    const s = stateRef.current;
    s.frame++;
    const cells = s.cells;

    // 1. Divide
    const toDiv = cells.filter(c => (s.frame - c.cycle.birthFrame) >= c.cycle.duration);
    if (toDiv.length && cells.length < MAX_CELLS) {
      const newCells = [];
      const divIds   = new Set(toDiv.map(c => c.id));
      for (const parent of toDiv) {
        if (cells.length + newCells.length >= MAX_CELLS) break;
        const angle = Math.random() * Math.PI * 2;
        const dx = Math.cos(angle) * BIRTH_SEP;
        const dy = Math.sin(angle) * BIRTH_SEP;
        newCells.push(makeCell(parent.x + dx, parent.y + dy, parent.generation + 1, s.frame));
        newCells.push(makeCell(parent.x - dx, parent.y - dy, parent.generation + 1, s.frame));
      }
      s.cells = cells.filter(c => !divIds.has(c.id)).concat(newCells);
    }

    // 2. Forces
    const { fx, fy } = computeForces(s.cells);

    // 3. Integrate + boundary
    for (let i = 0; i < s.cells.length; i++) {
      const c = s.cells[i];
      c.x += (DT / DAMPING) * fx[i];
      c.y += (DT / DAMPING) * fy[i];
      // Soft wall
      const pad = CELL_RADIUS;
      if (c.x < pad)            c.x += 0.5 * (pad - c.x);
      if (c.x > s.width - pad)  c.x -= 0.5 * (c.x - (s.width - pad));
      if (c.y < pad)            c.y += 0.5 * (pad - c.y);
      if (c.y > s.height - pad) c.y -= 0.5 * (c.y - (s.height - pad));

      // Scale-in animation
      if (c.scale < 1) c.scale = Math.min(1, c.scale + 0.04);
      // Flash timer
      if (c.flashTimer > 0) c.flashTimer--;
    }
  }, []);

  // ── Draw ──────────────────────────────────────────────────────────────────
  const draw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx  = canvas.getContext("2d");
    const s    = stateRef.current;
    const W    = s.width, H = s.height;

    // Background
    ctx.fillStyle = "#050810";
    ctx.fillRect(0, 0, W, H);

    // Subtle grid
    ctx.strokeStyle = "rgba(255,255,255,0.03)";
    ctx.lineWidth = 1;
    for (let x = 0; x < W; x += 40) { ctx.beginPath(); ctx.moveTo(x, 0); ctx.lineTo(x, H); ctx.stroke(); }
    for (let y = 0; y < H; y += 40) { ctx.beginPath(); ctx.moveTo(0, y); ctx.lineTo(W, y); ctx.stroke(); }

    // Spring connections (draw first, behind cells)
    ctx.lineWidth = 0.6;
    const cells = s.cells;
    for (let i = 0; i < cells.length; i++) {
      for (let j = i + 1; j < cells.length; j++) {
        const dx = cells[j].x - cells[i].x;
        const dy = cells[j].y - cells[i].y;
        const d2 = dx * dx + dy * dy;
        if (d2 > CUTOFF * CUTOFF) continue;
        const alpha = 0.12 * (1 - Math.sqrt(d2) / CUTOFF);
        ctx.strokeStyle = `rgba(180,220,255,${alpha})`;
        ctx.beginPath();
        ctx.moveTo(cells[i].x, cells[i].y);
        ctx.lineTo(cells[j].x, cells[j].y);
        ctx.stroke();
      }
    }

    // Cells
    for (const c of cells) {
      const color  = GEN_COLORS[c.generation % GEN_COLORS.length];
      const r      = CELL_RADIUS * c.scale;
      const age    = s.frame - c.cycle.birthFrame;
      const phase  = Math.min(1, age / c.cycle.duration); // 0→1 through cycle

      // Outer glow
      const glow = ctx.createRadialGradient(c.x, c.y, r * 0.2, c.x, c.y, r * 2.2);
      glow.addColorStop(0, color + "44");
      glow.addColorStop(1, "transparent");
      ctx.fillStyle = glow;
      ctx.beginPath();
      ctx.arc(c.x, c.y, r * 2.2, 0, Math.PI * 2);
      ctx.fill();

      // Cell body gradient
      const grad = ctx.createRadialGradient(c.x - r * 0.3, c.y - r * 0.3, r * 0.1, c.x, c.y, r);
      grad.addColorStop(0, color + "ff");
      grad.addColorStop(0.6, color + "cc");
      grad.addColorStop(1, color + "55");
      ctx.fillStyle = grad;
      ctx.beginPath();
      ctx.arc(c.x, c.y, r, 0, Math.PI * 2);
      ctx.fill();

      // Cell cycle ring — shows how close to division
      const arcLen = phase * Math.PI * 2;
      ctx.strokeStyle = phase > 0.85 ? "#ffffff" : color + "99";
      ctx.lineWidth   = phase > 0.85 ? 2.5 : 1.2;
      ctx.beginPath();
      ctx.arc(c.x, c.y, r + 2, -Math.PI / 2, -Math.PI / 2 + arcLen);
      ctx.stroke();

      // Nucleus
      ctx.fillStyle = "rgba(255,255,255,0.25)";
      ctx.beginPath();
      ctx.arc(c.x, c.y, r * 0.28, 0, Math.PI * 2);
      ctx.fill();
    }

    // HUD
    const maxG = cells.reduce((m, c) => Math.max(m, c.generation), 0);
    ctx.fillStyle = "rgba(0,0,0,0.55)";
    ctx.roundRect(12, 12, 180, 80, 8);
    ctx.fill();
    ctx.font = "bold 11px 'Courier New', monospace";
    ctx.fillStyle = "#4cc9f0";
    ctx.fillText(`CELLS     ${cells.length}`, 24, 34);
    ctx.fillStyle = "#69ff47";
    ctx.fillText(`FRAME     ${s.frame}`, 24, 52);
    ctx.fillStyle = "#ffd93d";
    ctx.fillText(`MAX GEN   ${maxG}`, 24, 70);
    ctx.fillStyle = "#ff6b6b";
    ctx.fillText(`TIME      ${(s.frame * DT / 60).toFixed(2)} h`, 24, 88);

    // Generation legend
    const legX = W - 130, legY = 14;
    ctx.fillStyle = "rgba(0,0,0,0.55)";
    ctx.roundRect(legX - 8, legY, 122, Math.min(maxG + 2, 10) * 16 + 10, 8);
    ctx.fill();
    for (let g = 0; g <= Math.min(maxG, 9); g++) {
      const col = GEN_COLORS[g % GEN_COLORS.length];
      ctx.fillStyle = col;
      ctx.beginPath();
      ctx.arc(legX + 5, legY + 16 + g * 16, 5, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillStyle = "#cdd9e5";
      ctx.font = "10px 'Courier New', monospace";
      ctx.fillText(`Gen ${g}`, legX + 16, legY + 20 + g * 16);
    }
  }, []);

  // ── Loop ──────────────────────────────────────────────────────────────────
  const loop = useCallback(() => {
    const s = stateRef.current;
    if (!s.running) return;
    for (let i = 0; i < s.speed; i++) step();
    draw();
    const atMax = s.cells.length >= MAX_CELLS;
    setUiState(u => ({
      ...u,
      frame: s.frame,
      count: s.cells.length,
      maxReached: atMax,
    }));
    if (atMax) {
      s.running = false;
      setUiState(u => ({ ...u, running: false }));
      return;
    }
    rafRef.current = requestAnimationFrame(loop);
  }, [step, draw]);

  const startStop = useCallback(() => {
    const s = stateRef.current;
    s.running = !s.running;
    setUiState(u => ({ ...u, running: s.running }));
    if (s.running) rafRef.current = requestAnimationFrame(loop);
  }, [loop]);

  const setSpeed = useCallback((v) => {
    stateRef.current.speed = v;
    setUiState(u => ({ ...u, speed: v }));
  }, []);

  const handleReset = useCallback(() => {
    cancelAnimationFrame(rafRef.current);
    stateRef.current.running = false;
    reset();
    draw();
    setUiState(u => ({ ...u, running: false, frame: 0, count: 1, maxReached: false }));
  }, [reset, draw]);

  // ── Resize canvas to container ────────────────────────────────────────────
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const container = canvas.parentElement;
    const size = Math.min(container.clientWidth, container.clientHeight, 620);
    canvas.width  = size;
    canvas.height = size;
    stateRef.current.width  = size;
    stateRef.current.height = size;
    reset();
    draw();
  }, [reset, draw]);

  // Cleanup
  useEffect(() => () => cancelAnimationFrame(rafRef.current), []);

  const { running, count, frame, speed, maxReached } = uiState;

  return (
    <div style={{
      minHeight: "100vh",
      background: "#020408",
      display: "flex",
      flexDirection: "column",
      alignItems: "center",
      justifyContent: "center",
      fontFamily: "'Courier New', monospace",
      padding: "16px",
      boxSizing: "border-box",
    }}>
      {/* Title */}
      <div style={{ marginBottom: 14, textAlign: "center" }}>
        <div style={{
          fontSize: 13,
          letterSpacing: "0.25em",
          color: "#4cc9f0",
          textTransform: "uppercase",
          marginBottom: 4,
        }}>Chaste · Off-Lattice Simulation</div>
        <div style={{
          fontSize: 10,
          color: "#445566",
          letterSpacing: "0.15em",
        }}>Meineke Spring Force · Stochastic Cell Cycle</div>
      </div>

      {/* Canvas */}
      <div style={{
        position: "relative",
        border: "1px solid #1a2a3a",
        borderRadius: 4,
        boxShadow: "0 0 40px rgba(76,201,240,0.08)",
        width: "min(620px, 100vw - 32px)",
        height: "min(620px, 100vw - 32px)",
      }}>
        <canvas
          ref={canvasRef}
          style={{ display: "block", width: "100%", height: "100%", borderRadius: 4 }}
        />
        {maxReached && (
          <div style={{
            position: "absolute", inset: 0,
            display: "flex", alignItems: "center", justifyContent: "center",
            background: "rgba(2,4,8,0.7)", borderRadius: 4,
            flexDirection: "column", gap: 8,
          }}>
            <div style={{ color: "#ffd93d", fontSize: 14, letterSpacing: "0.2em" }}>
              MAX CELLS REACHED
            </div>
            <div style={{ color: "#445566", fontSize: 10 }}>{count} cells · {(frame * DT / 60).toFixed(2)} h simulated</div>
          </div>
        )}
      </div>

      {/* Controls */}
      <div style={{
        marginTop: 16,
        display: "flex",
        gap: 10,
        alignItems: "center",
        flexWrap: "wrap",
        justifyContent: "center",
      }}>
        <button
          onClick={startStop}
          style={{
            background: running ? "#1a0a0a" : "#0a1a0a",
            border: `1px solid ${running ? "#ff6b6b" : "#69ff47"}`,
            color: running ? "#ff6b6b" : "#69ff47",
            padding: "7px 22px",
            fontSize: 11,
            letterSpacing: "0.2em",
            cursor: "pointer",
            borderRadius: 3,
            textTransform: "uppercase",
          }}
        >
          {running ? "⏸ Pause" : "▶ Run"}
        </button>

        <button
          onClick={handleReset}
          style={{
            background: "#0a0a1a",
            border: "1px solid #4cc9f0",
            color: "#4cc9f0",
            padding: "7px 22px",
            fontSize: 11,
            letterSpacing: "0.2em",
            cursor: "pointer",
            borderRadius: 3,
            textTransform: "uppercase",
          }}
        >
          ↺ Reset
        </button>

        <div style={{ display: "flex", gap: 4, alignItems: "center" }}>
          <span style={{ color: "#445566", fontSize: 10, letterSpacing: "0.15em" }}>SPEED</span>
          {[1, 2, 4, 8].map(v => (
            <button
              key={v}
              onClick={() => setSpeed(v)}
              style={{
                background: speed === v ? "#0d1a2a" : "transparent",
                border: `1px solid ${speed === v ? "#4cc9f0" : "#1a2a3a"}`,
                color: speed === v ? "#4cc9f0" : "#445566",
                padding: "5px 10px",
                fontSize: 10,
                cursor: "pointer",
                borderRadius: 3,
                letterSpacing: "0.1em",
              }}
            >
              {v}×
            </button>
          ))}
        </div>
      </div>

      {/* Footer note */}
      <div style={{
        marginTop: 12,
        fontSize: 9,
        color: "#2a3a4a",
        letterSpacing: "0.12em",
        textAlign: "center",
      }}>
        Ring = cell-cycle progress · White = near division · Max {MAX_CELLS} cells
      </div>
    </div>
  );
}
