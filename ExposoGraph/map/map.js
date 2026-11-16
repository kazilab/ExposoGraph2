// ═══════════════════════════════════════════════════════════════════════
// Knowledge Graph Visualization — D3.js Force-Directed Layout
// ═══════════════════════════════════════════════════════════════════════

// Streamlit custom component (st.components.v2), registered in
// ui_map_viewer.py. Python passes {graph, style} as `data`; every time that
// data changes, Streamlit calls the default export again on the same mounted
// element, and the map is rebuilt from a pristine copy of map.html's markup.

const FONTS_URL =
  "https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700"
  + "&family=JetBrains+Mono:wght@400;500&display=swap";

// Fonts must be registered on the page (not inside the component's shadow
// root) to be usable by it.
function ensureFonts() {
  if (document.querySelector(`link[href="${FONTS_URL}"]`)) return;
  const link = document.createElement("link");
  link.rel = "stylesheet";
  link.href = FONTS_URL;
  document.head.appendChild(link);
}

const instances = new WeakMap();

export default function (component) {
  const { parentElement, data } = component;
  ensureFonts();
  let inst = instances.get(parentElement);
  if (!inst) {
    inst = {
      template: parentElement.querySelector(".app").cloneNode(true),
      teardown: null,
      // Kept across re-renders of this mounted component (see renderMap).
      view: { layout: {}, zoom: null },
    };
    instances.set(parentElement, inst);
  }
  if (inst.teardown) inst.teardown();
  const root = inst.template.cloneNode(true);
  parentElement.querySelector(".app").replaceWith(root);
  inst.teardown = renderMap(root, data.graph, data.style, inst.view);
  return () => {
    if (inst.teardown) inst.teardown();
    inst.teardown = null;
  };
}

function renderMap(root, GRAPH_DATA, MAP_STYLE, view) {
  const byId = id => root.querySelector(`#${id}`);

  // Node/edge colors, legend labels, and sizes come from Python
  // (_app_shared.py) via MAP_STYLE, passed in alongside GRAPH_DATA.
  const NODE_COLORS = MAP_STYLE.node_colors;
  const NODE_LABELS = MAP_STYLE.node_labels;
  const EDGE_COLORS = MAP_STYLE.edge_colors;
  const EDGE_LABELS = MAP_STYLE.edge_labels;
  const NODE_SIZES = MAP_STYLE.node_sizes;

  // ── Setup ────────────────────────────────────────────────────────────
  const container = byId("graph-container");
  const svg = d3.select(byId("graph-svg"));
  const width = container.clientWidth;
  const height = container.clientHeight;

  // Deep clone data
  const nodes = GRAPH_DATA.nodes.map(d => ({...d}));
  const edges = GRAPH_DATA.edges.map(d => ({...d}));
  const DEFAULT_EDGE_OPACITY = d => d.type === "PATHWAY" ? 0.25 : 0.45;
  const DEFAULT_EDGE_WIDTH = d => d.type === "PATHWAY" ? 0.8 : 1.2;
  // Nodes/edges the Python side marked as below the active tissue threshold
  // (GraphEngine.dim_by_tissue_threshold) render at this reduced opacity
  // instead of being hidden, per the "gray out, don't remove" filter semantics.
  const DIMMED_OPACITY = 0.12;

  // ── Header Stats (dynamic, reflect current selection) ──────────────
  const statNodeCount = byId("node-count");
  const statEdgeCount = byId("edge-count");
  const statNodeTypeCount = byId("node-type-count");
  const statEdgeTypeCount = byId("edge-type-count");

  /**
   * Update the four header counters.
   *   visibleNodeIds:    Set<string> of node ids currently "active"; null = all.
   *   visibleEdgeIndices: Set<number> of edge array indices active; null = all.
   * Type counts reflect the distinct node/edge types present in the visible set.
   */
  function updateHeaderStats(visibleNodeIds = null, visibleEdgeIndices = null) {
    const activeNodes = visibleNodeIds
      ? nodes.filter(n => visibleNodeIds.has(n.id))
      : nodes;
    const activeEdges = visibleEdgeIndices
      ? edges.filter((_, i) => visibleEdgeIndices.has(i))
      : edges;
    const nodeTypes = new Set(activeNodes.map(n => n.type));
    const edgeTypes = new Set(activeEdges.map(e => e.type));
    statNodeCount.textContent = activeNodes.length;
    statEdgeCount.textContent = activeEdges.length;
    statNodeTypeCount.textContent = nodeTypes.size;
    statEdgeTypeCount.textContent = edgeTypes.size;
  }

  updateHeaderStats();

  // ── Build Legends ────────────────────────────────────────────────────
  /**
   * Build a tiny SVG glyph that mirrors the shape used for this node type on the
   * force-directed graph (diamond = Carcinogen, hexagon = DNA_Adduct,
   * rounded rect = Pathway/Tissue, circle = Enzyme/Gene/Metabolite).
   */
  function buildNodeLegendGlyph(type, color) {
    const size = 18;
    const c = size / 2;
    const r = 5;
    if (type === "Carcinogen") {
      const pts = `${c},${c - r * 1.3} ${c + r},${c} ${c},${c + r * 1.3} ${c - r},${c}`;
      return `<svg class="legend-shape" viewBox="0 0 ${size} ${size}">`
        + `<polygon points="${pts}" fill="${color}" fill-opacity="0.85" `
        + `stroke="${color}" stroke-width="1.5" stroke-opacity="0.5"/></svg>`;
    }
    if (type === "DNA_Adduct") {
      const hex = d3.range(6).map(i => {
        const a = (Math.PI / 3) * i - Math.PI / 6;
        return `${c + r * Math.cos(a)},${c + r * Math.sin(a)}`;
      }).join(" ");
      return `<svg class="legend-shape" viewBox="0 0 ${size} ${size}">`
        + `<polygon points="${hex}" fill="${color}" fill-opacity="0.85" `
        + `stroke="${color}" stroke-width="1.5" stroke-opacity="0.5"/></svg>`;
    }
    if (type === "Pathway") {
      const w = r * 2.4, h = r * 1.4;
      return `<svg class="legend-shape" viewBox="0 0 ${size} ${size}">`
        + `<rect x="${c - w / 2}" y="${c - h / 2}" width="${w}" height="${h}" rx="3" `
        + `fill="${color}" fill-opacity="0.5" `
        + `stroke="${color}" stroke-width="1" stroke-opacity="0.4"/></svg>`;
    }
    if (type === "Tissue") {
      const w = r * 2.2, h = r * 1.3;
      return `<svg class="legend-shape" viewBox="0 0 ${size} ${size}">`
        + `<rect x="${c - w / 2}" y="${c - h / 2}" width="${w}" height="${h}" rx="6" `
        + `fill="${color}" fill-opacity="0.72" `
        + `stroke="${color}" stroke-width="1.3" stroke-opacity="0.45"/></svg>`;
    }
    // Enzyme / Gene / Metabolite → circle, with the same fill-opacity ladder
    const fillOp = type === "Enzyme" ? 0.85 : type === "Gene" ? 0.78 : 0.7;
    return `<svg class="legend-shape" viewBox="0 0 ${size} ${size}">`
      + `<circle cx="${c}" cy="${c}" r="${r}" fill="${color}" fill-opacity="${fillOp}" `
      + `stroke="${color}" stroke-width="1.5" stroke-opacity="0.4"/></svg>`;
  }

  const nodeLegend = byId("node-legend");
  Object.entries(NODE_COLORS).forEach(([type, color]) => {
    const count = nodes.filter(n => n.type === type).length;
    if (count === 0) return;
    const item = document.createElement("div");
    item.className = "legend-item";
    item.dataset.type = type;
    item.innerHTML = `${buildNodeLegendGlyph(type, color)}
      <span class="legend-label">${NODE_LABELS[type]}</span>
      <span class="legend-count">${count}</span>`;
    nodeLegend.appendChild(item);
  });

  const edgeLegend = byId("edge-legend");
  Object.entries(EDGE_COLORS).forEach(([type, color]) => {
    const count = edges.filter(e => e.type === type).length;
    if (count === 0) return;
    const item = document.createElement("div");
    item.className = "edge-legend-item";
    item.innerHTML = `<span class="edge-line" style="background:${color}"></span>
      <span class="legend-label">${EDGE_LABELS[type]}</span>
      <span class="legend-count">${count}</span>`;
    edgeLegend.appendChild(item);
  });

  // Node-type, carcinogen-group, and tissue-threshold filtering all happen in
  // Python now (see GraphEngine.map_viewer_subgraph/dim_by_tissue_threshold and
  // ui_map_viewer.py) -- this document only ever receives the already-filtered
  // GRAPH_DATA it's rendering, plus a per-node/edge `_dimmed` flag for the
  // tissue-threshold overlay (applied below in resetGraphVisibility()).
  const headerSubtitle = byId("header-subtitle");
  const tissueSourceNote = "GTEx v8 nTPM via Human Protein Atlas v25; weights are per-gene max-normalized across the bundled tissue panel.";
  const tissueCaveatNote = "Colon averages sigmoid/transverse; Breast uses mammary tissue; Kidney uses cortex; Esophagus uses mucosa. GSTT1 has no bundled HPA/GTEx nTPM value.";

  function updateBundleSubtitle() {
    const heavyMetalDetected = nodes.some(
      n => n.group === "Heavy Metals" || ["As", "Cr_V", "Ni2+", "Cd", "Be"].includes(n.id),
    );
    const subtitleBits = [
      `${nodes.length} nodes / ${edges.length} edges`
    ];
    headerSubtitle.textContent = `Bundled reference graph · ${subtitleBits.join(" · ")}`;
  }

  updateBundleSubtitle();

  // ── D3 Force Simulation ─────────────────────────────────────────────
  const g = svg.append("g");

  // Zoom
  const zoom = d3.zoom()
    .scaleExtent([0.15, 5])
    .on("zoom", (e) => {
      g.attr("transform", e.transform);
      const { x, y, k } = e.transform;
      view.zoom = [x, y, k];
    });
  svg.call(zoom);

  // Arrow markers
  const defs = svg.append("defs");
  Object.entries(EDGE_COLORS).forEach(([type, color]) => {
    defs.append("marker")
      .attr("id", `arrow-${type}`)
      .attr("viewBox", "0 -5 10 10")
      .attr("refX", 20)
      .attr("refY", 0)
      .attr("markerWidth", 5)
      .attr("markerHeight", 5)
      .attr("orient", "auto")
      .append("path")
      .attr("d", "M0,-4L10,0L0,4")
      .attr("fill", color)
      .attr("opacity", 0.7);
  });

  // Links
  const link = g.append("g")
    .selectAll("line")
    .data(edges)
    .join("line")
    .attr("stroke", d => EDGE_COLORS[d.type] || "#555")
    .attr("stroke-width", d => d.type === "PATHWAY" ? 0.8 : 1.2)
    .attr("stroke-opacity", d => d.type === "PATHWAY" ? 0.25 : 0.45)
    .attr("stroke-dasharray", d => d.type === "SCHEMATIC_LINKAGE" ? "4,3" : null)
    .attr("marker-end", d => d.type !== "PATHWAY" ? `url(#arrow-${d.type})` : null);

  // Nodes
  const node = g.append("g")
    .selectAll("g")
    .data(nodes)
    .join("g")
    .call(d3.drag()
      .on("start", dragstarted)
      .on("drag", dragged)
      .on("end", dragended));

  // Node shapes
  node.each(function(d) {
    const el = d3.select(this);
    const r = NODE_SIZES[d.type] || 7;
    const color = NODE_COLORS[d.type] || "#888";

    if (d.type === "Carcinogen") {
      // Diamond shape for carcinogens
      el.append("polygon")
        .attr("points", `0,${-r*1.3} ${r},0 0,${r*1.3} ${-r},0`)
        .attr("fill", color)
        .attr("fill-opacity", 0.85)
        .attr("stroke", color)
        .attr("stroke-width", 1.5)
        .attr("stroke-opacity", 0.5);
    } else if (d.type === "DNA_Adduct") {
      // Hexagon for DNA adducts
      const hex = d3.range(6).map(i => {
        const a = (Math.PI / 3) * i - Math.PI / 6;
        return `${r * Math.cos(a)},${r * Math.sin(a)}`;
      }).join(" ");
      el.append("polygon")
        .attr("points", hex)
        .attr("fill", color)
        .attr("fill-opacity", 0.85)
        .attr("stroke", color)
        .attr("stroke-width", 1.5)
        .attr("stroke-opacity", 0.5);
    } else if (d.type === "Pathway") {
      // Rounded rect for pathways
      el.append("rect")
        .attr("x", -r * 1.2)
        .attr("y", -r * 0.7)
        .attr("width", r * 2.4)
        .attr("height", r * 1.4)
        .attr("rx", 3)
        .attr("fill", color)
        .attr("fill-opacity", 0.5)
        .attr("stroke", color)
        .attr("stroke-width", 1)
        .attr("stroke-opacity", 0.4);
    } else if (d.type === "Tissue") {
      // Rounded rect for tissues
      el.append("rect")
        .attr("x", -r * 1.1)
        .attr("y", -r * 0.65)
        .attr("width", r * 2.2)
        .attr("height", r * 1.3)
        .attr("rx", 6)
        .attr("fill", color)
        .attr("fill-opacity", 0.72)
        .attr("stroke", color)
        .attr("stroke-width", 1.3)
        .attr("stroke-opacity", 0.45);
    } else {
      // Circle for enzymes, genes, and metabolites
      el.append("circle")
        .attr("r", r)
        .attr("fill", color)
        .attr("fill-opacity", d.type === "Enzyme" ? 0.85 : d.type === "Gene" ? 0.78 : 0.7)
        .attr("stroke", color)
        .attr("stroke-width", 1.5)
        .attr("stroke-opacity", 0.4);
    }
  });

  // Node labels
  node.append("text")
    .text(d => d.label)
    .attr("font-size", d => d.type === "Carcinogen" ? "10px" : d.type === "Pathway" ? "7px" : "8px")
    .attr("font-weight", d => d.type === "Carcinogen" ? "600" : "400")
    .attr("font-family", "Inter, sans-serif")
    .attr("fill", "#2a2d33")
    .attr("text-anchor", "middle")
    .attr("dy", d => (NODE_SIZES[d.type] || 7) + 12)
    .attr("paint-order", "stroke")
    .attr("stroke", "var(--bg)")
    .attr("stroke-width", 3);

  // Apply the initial opacity baseline (full for ordinary nodes/edges, dimmed
  // for anything the Python side flagged via _dimmed) before the simulation
  // starts, rather than waiting for a search/hover interaction to trigger it.
  resetGraphVisibility();

  // ── Layout persistence ───────────────────────────────────────────────
  // Every Streamlit widget change rebuilds the map with new data, which would
  // otherwise restart the force layout from scratch. Node positions are kept in
  // the component's `view` state (keyed by node id) so a re-render restores the
  // previous layout instead of re-running it. Positions for nodes that are not
  // in the current subgraph are kept, so toggling them back restores them too.
  let tickCount = 0;

  function loadLayout() {
    return view.layout;
  }

  function saveLayout() {
    nodes.forEach(n => { view.layout[n.id] = [n.x, n.y]; });
  }

  const savedLayout = loadLayout();
  const restoredIds = new Set();
  nodes.forEach(n => {
    const p = savedLayout[n.id];
    if (p) { n.x = p[0]; n.y = p[1]; restoredIds.add(n.id); }
  });

  // ── Simulation ───────────────────────────────────────────────────────
  const simulation = d3.forceSimulation(nodes)
    .force("link", d3.forceLink(edges).id(d => d.id).distance(d => {
      if (d.type === "PATHWAY") return 160;
      if (d.type === "FORMS_ADDUCT") return 80;
      return 100;
    }).strength(d => d.type === "PATHWAY" ? 0.1 : 0.4))
    .force("charge", d3.forceManyBody().strength(d => {
      if (d.type === "Carcinogen") return -400;
      if (d.type === "Pathway") return -200;
      return -250;
    }))
    .force("center", d3.forceCenter(width / 2, height / 2))
    .force("collision", d3.forceCollide().radius(d => (NODE_SIZES[d.type] || 7) + 15))
    .force("x", d3.forceX(width / 2).strength(0.03))
    .force("y", d3.forceY(height / 2).strength(0.03))
    .alphaDecay(0.02)
    .on("tick", ticked)
    .on("end", saveLayout);

  if (restoredIds.size > 0) {
    // Re-render with a known layout: don't animate. Nodes new to this subgraph
    // are settled off-screen against the restored nodes (held fixed), then the
    // simulation stays stopped until the user drags a node.
    simulation.stop();
    const pending = nodes.filter(n => !restoredIds.has(n.id));
    if (pending.length > 0) {
      nodes.forEach(n => { if (restoredIds.has(n.id)) { n.fx = n.x; n.fy = n.y; } });
      simulation.alpha(1);
      for (let i = 0; i < 300; i++) simulation.tick();
      nodes.forEach(n => { n.fx = null; n.fy = null; });
      saveLayout();
    }
    nodes.forEach(n => { n.vx = 0; n.vy = 0; });
    ticked();
  }

  function ticked() {
    // Checkpoint the layout while it is still animating, so a widget change
    // before the first layout finishes still has positions to restore.
    if (++tickCount % 25 === 0) saveLayout();
    link
      .attr("x1", d => d.source.x)
      .attr("y1", d => d.source.y)
      .attr("x2", d => d.target.x)
      .attr("y2", d => d.target.y);

    node.attr("transform", d => `translate(${d.x},${d.y})`);
  }

  function resetGraphVisibility() {
    node.style("opacity", d => d._dimmed ? DIMMED_OPACITY : 1);
    link
      .style("stroke-opacity", d => d._dimmed ? DIMMED_OPACITY : DEFAULT_EDGE_OPACITY(d))
      .style("stroke-width", d => DEFAULT_EDGE_WIDTH(d));
    updateHeaderStats();
  }

  function normalizeString(value) {
    return String(value || "").toLowerCase();
  }

  function sortedTissueEntries(weights) {
    return Object.entries(weights || {}).sort((a, b) => b[1] - a[1]);
  }

  function formatTissueWeights(weights, preferredTissue = null, limit = 4) {
    let entries = sortedTissueEntries(weights);
    if (preferredTissue) {
      const preferred = entries.find(([name]) => name === preferredTissue);
      const others = entries.filter(([name]) => name !== preferredTissue);
      entries = preferred ? [preferred, ...others] : entries;
    }
    return entries
      .slice(0, limit)
      .map(([name, value]) => `${name} ${Number(value).toFixed(2)}`)
      .join(" · ");
  }

  // The only remaining client-side "view" is free-text search; everything else
  // (node type, carcinogen group, tissue threshold) is decided in Python before
  // this document is rendered. Re-applying the current search box value here
  // restores either the active search highlight or the tissue-dimmed baseline
  // (applySearch("") falls through to resetGraphVisibility()).
  function restoreView() {
    applySearch(searchInput.value.toLowerCase().trim());
  }

  // ── Drag ─────────────────────────────────────────────────────────────
  // d3.drag starts a gesture on every mousedown, including a plain click on a
  // node. The layout is only reheated once the pointer has actually moved the
  // node, so clicking a node (e.g. to open its details) never moves the graph.
  const DRAG_THRESHOLD_PX = 3;
  let dragMoved = false;
  function dragstarted() {
    dragMoved = false;
  }
  function dragged(event, d) {
    if (!dragMoved) {
      if (Math.hypot(event.x - d.x, event.y - d.y) < DRAG_THRESHOLD_PX) return;
      dragMoved = true;
      simulation.alphaTarget(0.3).restart();
    }
    d.fx = event.x;
    d.fy = event.y;
  }
  function dragended(event, d) {
    if (!dragMoved) return;
    simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
    saveLayout();
  }

  // ── Tooltip ──────────────────────────────────────────────────────────
  const tooltip = byId("tooltip");

  node.on("mouseover", (event, d) => {
    const color = NODE_COLORS[d.type];
    let typeExtra = "";
    if (d.type === "Enzyme") typeExtra = ` · Phase ${d.phase} · ${d.role}`;
    if (d.type === "Carcinogen") {
      const carcinogenBits = [
        d.group ? d.group.replace("_", " ") : null,
        d.iarc || null,
      ].filter(Boolean);
      if (carcinogenBits.length) typeExtra = ` · ${carcinogenBits.join(" · ")}`;
    }
    if (d.type === "Metabolite") typeExtra = ` · Reactivity: ${d.reactivity}`;
    if (d.type === "Gene") {
      const geneBits = [d.variant, d.phenotype, d.tissue].filter(Boolean);
      if (geneBits.length) typeExtra = ` · ${geneBits.join(" · ")}`;
    }
    if (d.type === "Tissue" && d.tissue) typeExtra = ` · ${d.tissue}`;
    const tissueWeightsHtml = d.tissue_weights
      ? `<div class="tooltip-desc"><b>GTEx/HPA nTPM weights:</b> ${formatTissueWeights(d.tissue_weights, null, 3)}</div>`
      : "";

    tooltip.innerHTML = `
      <div class="tooltip-title">${d.label}</div>
      <div class="tooltip-type" style="color:${color}">${d.type.replace("_"," ")}${typeExtra}</div>
      <div class="tooltip-desc">${d.detail}</div>
      ${tissueWeightsHtml}
    `;
    tooltip.classList.add("visible");

    // Highlight connected
    const connectedIds = new Set();
    connectedIds.add(d.id);
    edges.forEach(e => {
      const sid = typeof e.source === "object" ? e.source.id : e.source;
      const tid = typeof e.target === "object" ? e.target.id : e.target;
      if (sid === d.id) connectedIds.add(tid);
      if (tid === d.id) connectedIds.add(sid);
    });

    node.style("opacity", n => connectedIds.has(n.id) ? 1 : 0.1);
    link.style("stroke-opacity", e => {
      const sid = typeof e.source === "object" ? e.source.id : e.source;
      const tid = typeof e.target === "object" ? e.target.id : e.target;
      return (sid === d.id || tid === d.id) ? 0.8 : 0.03;
    }).style("stroke-width", e => {
      const sid = typeof e.source === "object" ? e.source.id : e.source;
      const tid = typeof e.target === "object" ? e.target.id : e.target;
      return (sid === d.id || tid === d.id) ? 2.5 : 0.5;
    });
  })
  .on("mousemove", (event) => {
    const rect = container.getBoundingClientRect();
    tooltip.style.left = (event.clientX - rect.left + 15) + "px";
    tooltip.style.top = (event.clientY - rect.top - 10) + "px";
  })
  .on("mouseout", () => {
    tooltip.classList.remove("visible");
    restoreView();
  });

  // ── Click for Details ────────────────────────────────────────────────
  node.on("click", (event, d) => {
    event.stopPropagation();
    showDetail(d);
  });

  svg.on("click", () => {
    byId("detail-panel").classList.remove("visible");
  });

  function showDetail(d) {
    const panel = byId("detail-panel");
    const content = byId("detail-content");
    const color = NODE_COLORS[d.type];

    let typeLabel = d.type.replace("_", " ");
    if (d.type === "Enzyme") typeLabel += ` · Phase ${d.phase}`;
    if (d.type === "Gene" && d.variant) typeLabel += ` · ${d.variant}`;

    // Find connections
    const outgoing = [];
    const incoming = [];
    edges.forEach(e => {
      const sid = typeof e.source === "object" ? e.source.id : e.source;
      const tid = typeof e.target === "object" ? e.target.id : e.target;
      if (sid === d.id) {
        const tNode = nodes.find(n => n.id === tid);
        outgoing.push({ edge: e, node: tNode });
      }
      if (tid === d.id) {
        const sNode = nodes.find(n => n.id === sid);
        incoming.push({ edge: e, node: sNode });
      }
    });

    let connHtml = "";
    if (outgoing.length) {
      connHtml += `<h4>Outgoing (${outgoing.length})</h4>`;
      outgoing.forEach(c => {
        const ec = EDGE_COLORS[c.edge.type] || "#555";
        connHtml += `<div class="conn-item">
          <span class="conn-edge" style="background:${ec}22;color:${ec}">${c.edge.type}</span>
          → ${c.node?.label || "?"}${c.edge.label ? ` (${c.edge.label})` : ""}
        </div>`;
      });
    }
    if (incoming.length) {
      connHtml += `<h4 style="margin-top:8px">Incoming (${incoming.length})</h4>`;
      incoming.forEach(c => {
        const ec = EDGE_COLORS[c.edge.type] || "#555";
        connHtml += `<div class="conn-item">
          <span class="conn-edge" style="background:${ec}22;color:${ec}">${c.edge.type}</span>
          ← ${c.node?.label || "?"}${c.edge.label ? ` (${c.edge.label})` : ""}
        </div>`;
      });
    }

    const metadata = [
      d.group ? `<div><b>Group:</b> ${d.group.replaceAll("_", " ")}</div>` : "",
      d.iarc ? `<div><b>IARC:</b> ${d.iarc}</div>` : "",
      d.phase ? `<div><b>Phase:</b> ${d.phase}</div>` : "",
      d.role ? `<div><b>Role:</b> ${d.role}</div>` : "",
      d.reactivity ? `<div><b>Reactivity:</b> ${d.reactivity}</div>` : "",
      d.tissue ? `<div><b>Tissue:</b> ${d.tissue}</div>` : "",
      d.tissue_weights ? `<div><b>GTEx/HPA nTPM Weights:</b> ${formatTissueWeights(d.tissue_weights, null, 8)}</div>` : "",
      d.tissue_weights ? `<div><b>Tissue Data:</b> ${tissueSourceNote} ${tissueCaveatNote}</div>` : "",
      d.variant ? `<div><b>Variant:</b> ${d.variant}</div>` : "",
      d.phenotype ? `<div><b>Phenotype:</b> ${d.phenotype}</div>` : "",
      d.source_db ? `<div><b>Source DB:</b> ${d.source_db}</div>` : "",
    ].filter(Boolean).join("");

    content.innerHTML = `
      <div class="detail-name">${d.label}</div>
      <span class="detail-type" style="background:${color}22;color:${color}">${typeLabel}</span>
      <div class="detail-desc">${d.detail}</div>
      ${metadata ? `<div class="detail-desc">${metadata}</div>` : ""}
      <div class="detail-connections">${connHtml}</div>
    `;
    panel.classList.add("visible");
  }

  // ── Search ───────────────────────────────────────────────────────────
  function applySearch(q) {
    if (!q) {
      resetGraphVisibility();
      return;
    }

    const matchIds = new Set();
    nodes.forEach(n => {
      const tissueWeightText = n.tissue_weights
        ? Object.entries(n.tissue_weights).map(([name, value]) => `${name} ${value}`).join(" ")
        : "";
      if (n.label.toLowerCase().includes(q) || n.id.toLowerCase().includes(q) ||
          (n.detail && n.detail.toLowerCase().includes(q)) ||
          (n.group && n.group.toLowerCase().includes(q)) ||
          (n.tissue && n.tissue.toLowerCase().includes(q)) ||
          tissueWeightText.toLowerCase().includes(q) ||
          (n.variant && n.variant.toLowerCase().includes(q)) ||
          (n.phenotype && n.phenotype.toLowerCase().includes(q))) {
        matchIds.add(n.id);
      }
    });

    // Also highlight connected nodes
    const expandedIds = new Set(matchIds);
    edges.forEach(e => {
      const sid = typeof e.source === "object" ? e.source.id : e.source;
      const tid = typeof e.target === "object" ? e.target.id : e.target;
      if (matchIds.has(sid)) expandedIds.add(tid);
      if (matchIds.has(tid)) expandedIds.add(sid);
    });

    const visibleEdgeIdx = new Set();
    edges.forEach((e, i) => {
      const sid = typeof e.source === "object" ? e.source.id : e.source;
      const tid = typeof e.target === "object" ? e.target.id : e.target;
      if (expandedIds.has(sid) && expandedIds.has(tid)) visibleEdgeIdx.add(i);
    });

    node.style("opacity", n => expandedIds.has(n.id) ? 1 : 0.08);
    link.style("stroke-opacity", (e, i) => visibleEdgeIdx.has(i) ? 0.6 : 0.02);

    updateHeaderStats(expandedIds, visibleEdgeIdx);
  }

  const searchInput = byId("search-input");
  searchInput.addEventListener("input", function() {
    applySearch(this.value.toLowerCase().trim());
  });

  // Node-type highlighting, carcinogen-group filtering, and tissue-threshold
  // filtering used to live here as client-side JS (applyTypeHighlight/
  // highlightType, applyCarcinogenFilter/filterByCarcinogen, applyTissueFilter/
  // filterByTissue). They have moved to Python -- GraphEngine.map_viewer_subgraph
  // and GraphEngine.dim_by_tissue_threshold in engine.py -- and are now driven by
  // Streamlit widgets in ui_map_viewer.py, which recompute the subgraph and
  // re-render this document rather than manipulating it client-side.

  // ── Sidebar Resizer ──────────────────────────────────────────────────
  const resizeObserver = (function setupSidebarResizer() {
    const app = root;
    const resizer = byId("sidebar-resizer");
    if (!app || !resizer) return;

    const MIN_W = 240;
    const MAX_W = 640;
    let dragging = false;

    function setWidth(px) {
      const clamped = Math.min(MAX_W, Math.max(MIN_W, px));
      app.style.setProperty("--sidebar-w", clamped + "px");
    }

    resizer.addEventListener("pointerdown", (e) => {
      dragging = true;
      resizer.classList.add("dragging");
      app.classList.add("is-resizing");
      resizer.setPointerCapture(e.pointerId);
    });
    resizer.addEventListener("pointermove", (e) => {
      if (!dragging) return;
      setWidth(e.clientX - app.getBoundingClientRect().left);
    });
    const endDrag = (e) => {
      if (!dragging) return;
      dragging = false;
      resizer.classList.remove("dragging");
      app.classList.remove("is-resizing");
      try { resizer.releasePointerCapture(e.pointerId); } catch (_) {}
      rebalanceForces();
    };
    resizer.addEventListener("pointerup", endDrag);
    resizer.addEventListener("pointercancel", endDrag);

    // Double-click to reset to the responsive default
    resizer.addEventListener("dblclick", () => {
      app.style.removeProperty("--sidebar-w");
      rebalanceForces();
    });

    // Re-center the force simulation when the graph container resizes.
    const ro = new ResizeObserver(() => rebalanceForces());
    ro.observe(container);
    return ro;
  })();

  let lastLayoutSize = [width, height];
  function rebalanceForces() {
    const w = container.clientWidth;
    const h = container.clientHeight;
    if (!w || !h) return;
    if (w === lastLayoutSize[0] && h === lastLayoutSize[1]) return;
    lastLayoutSize = [w, h];
    simulation.force("center", d3.forceCenter(w / 2, h / 2));
    simulation.force("x", d3.forceX(w / 2).strength(0.03));
    simulation.force("y", d3.forceY(h / 2).strength(0.03));
    simulation.alpha(0.15).restart();
  }

  // ── Zoom Controls ────────────────────────────────────────────────────
  byId("zoom-in").addEventListener("click", () => {
    svg.transition().duration(300).call(zoom.scaleBy, 1.5);
  });
  byId("zoom-out").addEventListener("click", () => {
    svg.transition().duration(300).call(zoom.scaleBy, 0.67);
  });
  byId("zoom-reset").addEventListener("click", () => {
    svg.transition().duration(500).call(zoom.transform, d3.zoomIdentity);
  });

  // Initial fit
  const savedZoom = view.zoom;
  let initialFitTimer = null;
  if (savedZoom) {
    const [x, y, k] = savedZoom;
    svg.call(zoom.transform, d3.zoomIdentity.translate(x, y).scale(k));
  } else {
    initialFitTimer = setTimeout(() => {
      svg.transition().duration(800).call(
        zoom.transform,
        d3.zoomIdentity.translate(width * 0.05, height * 0.05).scale(0.85)
      );
    }, 2000);
  }

  // Stop everything this render started, so the next render (new data from
  // Python) or an unmount doesn't leave a running simulation or observers.
  return () => {
    simulation.stop();
    resizeObserver?.disconnect();
    clearTimeout(initialFitTimer);
    svg.interrupt();
  };
}
