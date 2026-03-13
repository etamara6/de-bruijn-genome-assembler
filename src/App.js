import { useState, useCallback, useRef, useEffect } from "react";

// ─── Constants ────────────────────────────────────────────────────────────────

const API_URL = process.env.REACT_APP_API_URL || "";

const DEMO_READS = [
  "ACGTTGCATGC",
  "TGCATGCAAC",
];

const DNA_COLORS = { A: "#4ade80", C: "#60a5fa", G: "#f59e0b", T: "#f87171", N: "#9ca3af" };

// ─── Utility helpers ──────────────────────────────────────────────────────────

function validateReads(reads) {
  const errors = [];
  reads.forEach((read, i) => {
    const trimmed = read.trim();
    if (!trimmed) return;
    if (!/^[ACGTNacgtn]+$/.test(trimmed)) {
      errors.push(`Read ${i + 1}: invalid characters (only A, C, G, T, N allowed)`);
    }
    if (trimmed.length < 2) {
      errors.push(`Read ${i + 1}: too short (minimum 2 bp)`);
    }
  });
  return errors;
}

function formatMs(ms) {
  if (ms === null || ms === undefined) return "—";
  return ms < 1 ? `${(ms * 1000).toFixed(0)} µs` : `${ms.toFixed(2)} ms`;
}

function colorDNA(seq) {
  return seq.split("").map((base, i) => (
    <span key={i} style={{ color: DNA_COLORS[base.toUpperCase()] ?? "#e2e8f0" }}>
      {base}
    </span>
  ));
}

// ─── Sub-components ───────────────────────────────────────────────────────────

function Helix({ className = "", style: extraStyle = {} }) {
  return (
    <svg
      className={className}
      style={extraStyle}
      viewBox="0 0 60 200"
      fill="none"
      xmlns="http://www.w3.org/2000/svg"
      aria-hidden="true"
    >
      {Array.from({ length: 10 }).map((_, i) => {
        const y = 10 + i * 18;
        const phase = (i / 10) * Math.PI * 2;
        const x1 = 10 + 20 * Math.sin(phase);
        const x2 = 50 - 20 * Math.sin(phase);
        return (
          <g key={i}>
            <circle cx={x1} cy={y} r="3.5" fill="#34d399" opacity="0.85" />
            <circle cx={x2} cy={y} r="3.5" fill="#60a5fa" opacity="0.85" />
            <line x1={x1} y1={y} x2={x2} y2={y} stroke="#94a3b8" strokeWidth="1" opacity="0.4" />
          </g>
        );
      })}
      <path
        d={Array.from({ length: 10 }).map((_, i) => {
          const y = 10 + i * 18;
          const x = 10 + 20 * Math.sin((i / 10) * Math.PI * 2);
          return `${i === 0 ? "M" : "L"} ${x} ${y}`;
        }).join(" ")}
        stroke="#34d399" strokeWidth="1.5" fill="none" opacity="0.5" strokeDasharray="3 2"
      />
      <path
        d={Array.from({ length: 10 }).map((_, i) => {
          const y = 10 + i * 18;
          const x = 50 - 20 * Math.sin((i / 10) * Math.PI * 2);
          return `${i === 0 ? "M" : "L"} ${x} ${y}`;
        }).join(" ")}
        stroke="#60a5fa" strokeWidth="1.5" fill="none" opacity="0.5" strokeDasharray="3 2"
      />
    </svg>
  );
}

function StatCard({ label, value, accent }) {
  return (
    <div style={{
      background: "rgba(15,23,42,0.6)",
      border: "1px solid rgba(148,163,184,0.12)",
      borderRadius: "10px",
      padding: "14px 18px",
      display: "flex",
      flexDirection: "column",
      gap: "4px",
      minWidth: "110px",
    }}>
      <span style={{ fontSize: "11px", color: "#64748b", fontFamily: "var(--font-mono)", letterSpacing: "0.08em", textTransform: "uppercase" }}>
        {label}
      </span>
      <span style={{ fontSize: "20px", fontWeight: "700", color: accent || "#e2e8f0", fontFamily: "var(--font-mono)" }}>
        {value}
      </span>
    </div>
  );
}

function ContigCard({ contig, index }) {
  const [copied, setCopied] = useState(false);

  const copy = useCallback(() => {
    navigator.clipboard.writeText(contig).then(() => {
      setCopied(true);
      setTimeout(() => setCopied(false), 1500);
    });
  }, [contig]);

  return (
    <div style={{
      background: "rgba(15,23,42,0.7)",
      border: "1px solid rgba(52,211,153,0.18)",
      borderRadius: "12px",
      padding: "18px 20px",
      display: "flex",
      flexDirection: "column",
      gap: "10px",
    }}>
      <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between" }}>
        <span style={{ fontSize: "11px", color: "#34d399", fontFamily: "var(--font-mono)", letterSpacing: "0.1em" }}>
          CONTIG_{String(index + 1).padStart(2, "0")}  ·  {contig.length} bp
        </span>
        <button
          onClick={copy}
          style={{
            background: copied ? "rgba(52,211,153,0.15)" : "transparent",
            border: "1px solid rgba(148,163,184,0.2)",
            borderRadius: "6px",
            color: copied ? "#34d399" : "#94a3b8",
            cursor: "pointer",
            fontSize: "11px",
            fontFamily: "var(--font-mono)",
            padding: "3px 10px",
            transition: "all 0.15s",
          }}
        >
          {copied ? "✓ copied" : "copy"}
        </button>
      </div>
      <div style={{
        fontFamily: "var(--font-mono)",
        fontSize: "13.5px",
        letterSpacing: "0.06em",
        lineHeight: "1.8",
        wordBreak: "break-all",
        padding: "10px 14px",
        background: "rgba(0,0,0,0.3)",
        borderRadius: "8px",
        borderLeft: "3px solid rgba(52,211,153,0.4)",
      }}>
        {colorDNA(contig)}
      </div>
    </div>
  );
}

function KmerPill({ kmer }) {
  return (
    <span style={{
      display: "inline-flex",
      fontFamily: "var(--font-mono)",
      fontSize: "11px",
      background: "rgba(96,165,250,0.08)",
      border: "1px solid rgba(96,165,250,0.2)",
      borderRadius: "5px",
      padding: "2px 8px",
      color: "#93c5fd",
      letterSpacing: "0.05em",
    }}>
      {kmer}
    </span>
  );
}

function GraphViz({ stats }) {
  if (!stats) return null;
  const { graph_nodes, graph_edges, total_kmers, n50 } = stats;
  const coverage = graph_edges > 0 ? ((graph_edges / total_kmers) * 100).toFixed(0) : 0;

  return (
    <div style={{
      background: "rgba(15,23,42,0.5)",
      border: "1px solid rgba(148,163,184,0.1)",
      borderRadius: "12px",
      padding: "18px 20px",
    }}>
      <div style={{ fontSize: "11px", color: "#64748b", fontFamily: "var(--font-mono)", letterSpacing: "0.1em", marginBottom: "14px" }}>
        GRAPH TOPOLOGY
      </div>
      <div style={{ display: "flex", flexWrap: "wrap", gap: "8px" }}>
        <StatCard label="nodes" value={graph_nodes} accent="#60a5fa" />
        <StatCard label="edges" value={graph_edges} accent="#a78bfa" />
        <StatCard label="k-mers" value={total_kmers} accent="#f59e0b" />
        <StatCard label="N50" value={n50 ? `${n50} bp` : "—"} accent="#34d399" />
        <StatCard label="coverage" value={`${coverage}%`} accent="#f472b6" />
      </div>
    </div>
  );
}

// ─── Main App ─────────────────────────────────────────────────────────────────

export default function App() {
  const [readsInput, setReadsInput] = useState(DEMO_READS.join("\n"));
  const [kValue, setKValue] = useState(4);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [result, setResult] = useState(null);
  const [validationErrors, setValidationErrors] = useState([]);
  const [apiStatus, setApiStatus] = useState("unknown");
  const resultRef = useRef(null);

  useEffect(() => {
    fetch(`${API_URL}/api/health`)
      .then((r) => (r.ok ? setApiStatus("ok") : setApiStatus("error")))
      .catch(() => setApiStatus("error"));
  }, []);

  useEffect(() => {
    const reads = readsInput.split("\n").filter((r) => r.trim());
    setValidationErrors(reads.length ? validateReads(reads) : []);
  }, [readsInput]);

  const handleAssemble = useCallback(async () => {
    setError(null);
    setResult(null);
    const reads = readsInput.split("\n").map((r) => r.trim()).filter(Boolean);
    if (!reads.length) { setError("Please enter at least one DNA read."); return; }
    const errs = validateReads(reads);
    if (errs.length) { setError(errs.join("\n")); return; }
    setLoading(true);
    try {
      const response = await fetch(`${API_URL}/api/assemble`, {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ reads, k: kValue }),
      });
      const data = await response.json();
      if (!response.ok) throw new Error(data.error || `Server error (${response.status})`);
      setResult(data);
      setTimeout(() => resultRef.current?.scrollIntoView({ behavior: "smooth", block: "start" }), 100);
    } catch (err) {
      setError(err.message || "Network error — is the Flask API running?");
    } finally {
      setLoading(false);
    }
  }, [readsInput, kValue]);

  const handleLoadDemo = useCallback(() => {
    setReadsInput(DEMO_READS.join("\n")); setKValue(4); setResult(null); setError(null);
  }, []);

  const handleClear = useCallback(() => {
    setReadsInput(""); setResult(null); setError(null);
  }, []);

  const reads = readsInput.split("\n").filter((r) => r.trim());
  const canSubmit = reads.length > 0 && validationErrors.length === 0 && !loading;

  return (
    <div style={{ minHeight: "100vh", background: "#030712", color: "#e2e8f0", fontFamily: "var(--font-body)", position: "relative", overflowX: "hidden" }}>
      <style>{`
        @import url('https://fonts.googleapis.com/css2?family=DM+Mono:ital,wght@0,300;0,400;0,500;1,300&family=Syne:wght@400;500;600;700;800&display=swap');
        :root { --font-body: 'Syne', sans-serif; --font-mono: 'DM Mono', monospace; }
        *, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
        body { background: #030712; }
        ::-webkit-scrollbar { width: 6px; }
        ::-webkit-scrollbar-track { background: #0f172a; }
        ::-webkit-scrollbar-thumb { background: #334155; border-radius: 3px; }
        textarea { resize: vertical; outline: none; }
        @keyframes pulse-dot { 0%,100% { opacity:0.6; transform:scale(0.95); } 50% { opacity:1; transform:scale(1.1); } }
        @keyframes slide-up { from { opacity:0; transform:translateY(20px); } to { opacity:1; transform:translateY(0); } }
        @keyframes spin { to { transform: rotate(360deg); } }
        .assemble-btn { position:relative; cursor:pointer; border:none; border-radius:10px; padding:13px 28px; font-family:var(--font-body); font-weight:700; font-size:15px; letter-spacing:0.05em; transition:transform 0.15s, box-shadow 0.15s, opacity 0.15s; }
        .assemble-btn:enabled:hover { transform:translateY(-2px); box-shadow:0 8px 30px rgba(52,211,153,0.25); }
        .assemble-btn:enabled:active { transform:translateY(0); }
        .assemble-btn:disabled { opacity:0.45; cursor:not-allowed; }
        .reads-textarea:focus { border-color:rgba(52,211,153,0.5)!important; box-shadow:0 0 0 3px rgba(52,211,153,0.08); }
        .result-section { animation: slide-up 0.4s ease both; }
        .k-input { outline:none; -moz-appearance:textfield; }
        .k-input::-webkit-inner-spin-button, .k-input::-webkit-outer-spin-button { -webkit-appearance:none; }
        .ghost-btn:hover { background:rgba(148,163,184,0.08)!important; }
        .demo-btn:hover { background:rgba(96,165,250,0.1)!important; }
        details summary::-webkit-details-marker { display:none; }
      `}</style>

      {/* Background grid */}
      <div style={{ position:"fixed", inset:0, backgroundImage:`linear-gradient(rgba(52,211,153,0.025) 1px, transparent 1px), linear-gradient(90deg, rgba(52,211,153,0.025) 1px, transparent 1px)`, backgroundSize:"44px 44px", pointerEvents:"none", zIndex:0 }} />
      <div style={{ position:"fixed", top:"-15%", left:"50%", transform:"translateX(-50%)", width:"800px", height:"400px", background:"radial-gradient(ellipse, rgba(52,211,153,0.035) 0%, transparent 68%)", pointerEvents:"none", zIndex:0 }} />

      <div style={{ position:"relative", zIndex:1, maxWidth:"860px", margin:"0 auto", padding:"0 24px 80px" }}>

        {/* Header */}
        <header style={{ paddingTop:"56px", paddingBottom:"48px", display:"flex", alignItems:"flex-start", gap:"20px" }}>
          <Helix style={{ width:"48px", flexShrink:0, marginTop:"4px" }} />
          <div>
            <div style={{ display:"flex", alignItems:"center", gap:"10px", marginBottom:"10px", flexWrap:"wrap" }}>
              <span style={{ fontFamily:"var(--font-mono)", fontSize:"10px", color:"#34d399", letterSpacing:"0.18em", textTransform:"uppercase", background:"rgba(52,211,153,0.08)", border:"1px solid rgba(52,211,153,0.2)", borderRadius:"4px", padding:"3px 8px" }}>
                De Bruijn · Eulerian Path · Hierholzer O(E)
              </span>
              <span style={{ display:"flex", alignItems:"center", gap:"5px", fontFamily:"var(--font-mono)", fontSize:"10px", color: apiStatus==="ok" ? "#34d399" : apiStatus==="error" ? "#f87171" : "#64748b", letterSpacing:"0.1em" }}>
                <span style={{ width:"6px", height:"6px", borderRadius:"50%", background: apiStatus==="ok"?"#34d399": apiStatus==="error"?"#f87171":"#64748b", display:"inline-block", animation: apiStatus==="ok"?"pulse-dot 2.5s ease infinite":"none" }} />
                {apiStatus==="ok" ? "API online" : apiStatus==="error" ? "API offline" : "checking…"}
              </span>
            </div>
            <h1 style={{ fontSize:"clamp(26px,5vw,42px)", fontWeight:"800", lineHeight:"1.1", letterSpacing:"-0.025em", color:"#f8fafc", marginBottom:"10px" }}>
              Genome Assembler
            </h1>
            <p style={{ color:"#64748b", fontSize:"14.5px", lineHeight:"1.65", maxWidth:"480px" }}>
              Reconstruct DNA sequences from short reads using De Bruijn graphs and Hierholzer's algorithm. Enter reads below, choose k, and assemble.
            </p>
          </div>
        </header>

        {/* Base legend */}
        <div style={{ display:"flex", gap:"8px", marginBottom:"28px", flexWrap:"wrap", alignItems:"center" }}>
          {Object.entries(DNA_COLORS).map(([base, color]) => (
            <span key={base} style={{ display:"inline-flex", alignItems:"center", gap:"5px", padding:"2px 9px", borderRadius:"5px", background:`${color}12`, fontFamily:"var(--font-mono)", fontSize:"12px", fontWeight:"500", color }}>
              <span style={{ width:"5px", height:"5px", borderRadius:"50%", background:color }} />{base}
            </span>
          ))}
          <span style={{ fontSize:"11px", color:"#334155", fontFamily:"var(--font-mono)" }}>base colour coding</span>
        </div>

        {/* Input panel */}
        <section style={{ background:"rgba(15,23,42,0.7)", border:"1px solid rgba(148,163,184,0.1)", borderRadius:"16px", padding:"28px", backdropFilter:"blur(12px)", marginBottom:"20px" }}>
          <div style={{ display:"flex", alignItems:"center", justifyContent:"space-between", marginBottom:"20px" }}>
            <div>
              <h2 style={{ fontSize:"16px", fontWeight:"700", color:"#f1f5f9", letterSpacing:"-0.01em" }}>Input Reads</h2>
              <span style={{ fontSize:"12px", color:"#475569", fontFamily:"var(--font-mono)" }}>one sequence per line · A C G T N</span>
            </div>
            <div style={{ display:"flex", gap:"8px" }}>
              <button className="demo-btn" onClick={handleLoadDemo} style={{ background:"transparent", border:"1px solid rgba(96,165,250,0.3)", borderRadius:"7px", color:"#60a5fa", cursor:"pointer", fontSize:"12px", fontFamily:"var(--font-mono)", padding:"6px 14px", transition:"background 0.15s" }}>
                demo
              </button>
              <button className="ghost-btn" onClick={handleClear} style={{ background:"transparent", border:"1px solid rgba(148,163,184,0.15)", borderRadius:"7px", color:"#64748b", cursor:"pointer", fontSize:"12px", fontFamily:"var(--font-mono)", padding:"6px 14px", transition:"background 0.15s" }}>
                clear
              </button>
            </div>
          </div>

          <textarea
            className="reads-textarea"
            value={readsInput}
            onChange={(e) => setReadsInput(e.target.value)}
            placeholder={"ACGTTGCATGC\nTGCATGCAAC\nACGTACGT"}
            rows={6}
            spellCheck={false}
            style={{ width:"100%", background:"rgba(0,0,0,0.35)", border:"1px solid rgba(148,163,184,0.12)", borderRadius:"10px", color:"#e2e8f0", fontFamily:"var(--font-mono)", fontSize:"14px", letterSpacing:"0.05em", lineHeight:"1.9", padding:"14px 16px", transition:"border-color 0.2s, box-shadow 0.2s" }}
          />

          {reads.length > 0 && (
            <div style={{ display:"flex", gap:"6px", flexWrap:"wrap", marginTop:"10px", alignItems:"center" }}>
              {reads.map((r, i) => {
                const hasErr = validationErrors.some(e => e.startsWith(`Read ${i + 1}:`));
                return (
                  <span key={i} style={{ fontFamily:"var(--font-mono)", fontSize:"11px", background: hasErr?"rgba(248,113,113,0.09)":"rgba(52,211,153,0.07)", border:`1px solid ${hasErr?"rgba(248,113,113,0.25)":"rgba(52,211,153,0.15)"}`, borderRadius:"5px", padding:"2px 8px", color: hasErr?"#fca5a5":"#6ee7b7" }}>
                    {r.trim().length} bp
                  </span>
                );
              })}
              <span style={{ fontSize:"11px", color:"#475569", fontFamily:"var(--font-mono)" }}>{reads.length} read{reads.length!==1?"s":""}</span>
            </div>
          )}

          {validationErrors.length > 0 && (
            <div style={{ marginTop:"12px", padding:"10px 14px", background:"rgba(248,113,113,0.07)", border:"1px solid rgba(248,113,113,0.2)", borderRadius:"8px" }}>
              {validationErrors.map((e, i) => (
                <div key={i} style={{ fontSize:"12px", color:"#fca5a5", fontFamily:"var(--font-mono)", lineHeight:"1.7" }}>⚠ {e}</div>
              ))}
            </div>
          )}

          <div style={{ display:"flex", alignItems:"center", gap:"16px", marginTop:"22px", flexWrap:"wrap" }}>
            {/* k stepper */}
            <div style={{ display:"flex", alignItems:"center", gap:"10px" }}>
              <label style={{ fontSize:"13px", color:"#94a3b8", fontFamily:"var(--font-mono)", letterSpacing:"0.05em" }}>k =</label>
              <div style={{ display:"flex", alignItems:"center", background:"rgba(0,0,0,0.35)", border:"1px solid rgba(148,163,184,0.15)", borderRadius:"8px", overflow:"hidden" }}>
                <button onClick={() => setKValue(v => Math.max(2, v-1))} style={{ background:"transparent", border:"none", color:"#94a3b8", cursor:"pointer", fontSize:"16px", padding:"6px 12px" }}>−</button>
                <input className="k-input" type="number" min="2" max="99" value={kValue} onChange={(e) => setKValue(Math.max(2, parseInt(e.target.value,10)||2))}
                  style={{ background:"transparent", border:"none", color:"#f59e0b", fontFamily:"var(--font-mono)", fontSize:"15px", fontWeight:"500", textAlign:"center", width:"42px", padding:"6px 0" }} />
                <button onClick={() => setKValue(v => Math.min(99, v+1))} style={{ background:"transparent", border:"none", color:"#94a3b8", cursor:"pointer", fontSize:"16px", padding:"6px 12px" }}>+</button>
              </div>
              <span style={{ fontSize:"11px", color:"#475569", fontFamily:"var(--font-mono)" }}>
                {kValue<=6?"testing": kValue<=15?"short reads":"real genome"}
              </span>
            </div>

            <button className="assemble-btn" onClick={handleAssemble} disabled={!canSubmit}
              style={{ background: canSubmit?"linear-gradient(135deg,#34d399 0%,#059669 100%)":"rgba(52,211,153,0.18)", color: canSubmit?"#022c22":"#34d399", marginLeft:"auto", minWidth:"150px" }}>
              {loading ? (
                <span style={{ display:"flex", alignItems:"center", gap:"8px", justifyContent:"center" }}>
                  <span style={{ width:"14px", height:"14px", border:"2px solid rgba(2,44,34,0.3)", borderTopColor:"#022c22", borderRadius:"50%", display:"inline-block", animation:"spin 0.7s linear infinite" }} />
                  assembling…
                </span>
              ) : "▶  Assemble"}
            </button>
          </div>
        </section>

        {/* Error */}
        {error && (
          <div style={{ background:"rgba(248,113,113,0.08)", border:"1px solid rgba(248,113,113,0.25)", borderRadius:"10px", padding:"14px 18px", marginBottom:"20px", fontFamily:"var(--font-mono)", fontSize:"13px", color:"#fca5a5", lineHeight:"1.6", animation:"slide-up 0.3s ease both" }}>
            <strong style={{ display:"block", marginBottom:"4px", color:"#f87171" }}>Error</strong>
            {error.split("\n").map((line, i) => <div key={i}>{line}</div>)}
          </div>
        )}

        {/* Results */}
        {result && (
          <div ref={resultRef} className="result-section">
            <div style={{ display:"flex", alignItems:"center", gap:"12px", marginBottom:"16px", flexWrap:"wrap" }}>
              <span style={{ fontFamily:"var(--font-mono)", fontSize:"10px", color:"#34d399", letterSpacing:"0.15em", textTransform:"uppercase" }}>Assembly complete</span>
              <span style={{ height:"1px", flex:1, background:"rgba(52,211,153,0.15)", minWidth:"20px" }} />
              <span style={{ fontFamily:"var(--font-mono)", fontSize:"11px", color:"#475569" }}>{formatMs(result.stats.elapsed_ms)}</span>
            </div>

            <GraphViz stats={result.stats} />

            <div style={{ marginTop:"16px" }}>
              <div style={{ display:"flex", alignItems:"center", justifyContent:"space-between", marginBottom:"14px" }}>
                <h3 style={{ fontSize:"15px", fontWeight:"700", color:"#f1f5f9" }}>Assembled Contigs</h3>
                <span style={{ fontFamily:"var(--font-mono)", fontSize:"11px", color:"#64748b", background:"rgba(148,163,184,0.07)", border:"1px solid rgba(148,163,184,0.1)", borderRadius:"5px", padding:"2px 8px" }}>
                  {result.contigs.length} contig{result.contigs.length!==1?"s":""} · longest-first
                </span>
              </div>

              {result.contigs.length === 0 ? (
                <div style={{ padding:"32px", textAlign:"center", background:"rgba(15,23,42,0.5)", border:"1px dashed rgba(148,163,184,0.15)", borderRadius:"12px", color:"#475569", fontFamily:"var(--font-mono)", fontSize:"13px" }}>
                  No contigs assembled — reads may not form a valid Eulerian path.<br />
                  <span style={{ fontSize:"11px", display:"block", marginTop:"6px" }}>Try different reads or a smaller k.</span>
                </div>
              ) : (
                <div style={{ display:"flex", flexDirection:"column", gap:"12px" }}>
                  {result.contigs.map((contig, i) => <ContigCard key={i} contig={contig} index={i} />)}
                </div>
              )}
            </div>

            {result.contigs.length > 0 && (
              <details style={{ marginTop:"16px" }}>
                <summary style={{ cursor:"pointer", fontFamily:"var(--font-mono)", fontSize:"11px", color:"#64748b", letterSpacing:"0.08em", userSelect:"none", padding:"10px 0", display:"flex", alignItems:"center", gap:"8px" }}>
                  <span style={{ fontSize:"10px" }}>▶</span>
                  K-MER PREVIEW  ·  k={kValue}  ·  first contig
                </summary>
                <div style={{ marginTop:"10px", padding:"14px 16px", background:"rgba(0,0,0,0.3)", border:"1px solid rgba(148,163,184,0.08)", borderRadius:"10px", display:"flex", flexWrap:"wrap", gap:"5px" }}>
                  {result.contigs[0].length >= kValue
                    ? Array.from({ length: result.contigs[0].length - kValue + 1 }, (_, i) => result.contigs[0].slice(i, i+kValue))
                        .map((kmer, i) => <KmerPill key={i} kmer={kmer} />)
                    : <span style={{ fontFamily:"var(--font-mono)", fontSize:"12px", color:"#475569" }}>Contig shorter than k.</span>
                  }
                </div>
              </details>
            )}
          </div>
        )}

        {/* Footer */}
        <footer style={{ marginTop:"72px", paddingTop:"24px", borderTop:"1px solid rgba(148,163,184,0.08)", display:"flex", alignItems:"center", justifyContent:"space-between", flexWrap:"wrap", gap:"12px" }}>
          <div style={{ display:"flex", alignItems:"center", gap:"10px" }}>
            <Helix style={{ width:"22px", opacity:0.4 }} />
            <span style={{ fontFamily:"var(--font-mono)", fontSize:"11px", color:"#334155" }}>
              Genome Assembler · De Bruijn Graph + Hierholzer O(E)
            </span>
          </div>
          <div style={{ display:"flex", gap:"16px" }}>
            {[["Python","#3b82f6"],["Flask","#94a3b8"],["React","#61dafb"]].map(([label,color])=>(
              <span key={label} style={{ fontFamily:"var(--font-mono)", fontSize:"10px", color, letterSpacing:"0.08em" }}>{label}</span>
            ))}
          </div>
        </footer>
      </div>
    </div>
  );
}