"use strict";

// The browser edits only the versioned request.  It does not calculate a case
// file or authorize a quantity: /api/preview and /api/compile always invoke the
// same Python compiler used by the command line.
let model;
let capabilities;
let materialLibrary;
let serverCapability;

const $ = (id) => document.getElementById(id);
const clone = (value) => JSON.parse(JSON.stringify(value));
const mm = (value) => ({value, unit:"mm"});
const identity = (z) => ({
  translation:{x:mm(0), y:mm(0), z:mm(z)},
  rotation_quaternion_local_to_world:{w:1, x:0, y:0, z:0}
});

async function api(path, body) {
  const options = body === undefined ? {} : {
    method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify(body)
  };
  const response = await fetch(path, options);
  const value = await response.json();
  if (!response.ok) throw new Error(value.error || `${response.status} ${response.statusText}`);
  return value;
}

function materialRecords() {
  const records = materialLibrary.materials.map(clone);
  for (const override of (model.material_overrides || [])) {
    if (override.extends) {
      const parentId=override.extends.split("@")[0];
      const parent=records.find((record)=>record.material_id===parentId);
      records.push({...clone(parent || {}),...clone(override),aliases:override.aliases || []});
    } else records.push(clone(override));
  }
  return records;
}

function materialsFor(role) {
  return materialRecords().filter((record) => (record.roles || []).includes(role));
}

function materialOptions(role, selected) {
  return materialsFor(role).map((record) => {
    const disabled = record.status === "template" ? " disabled" : "";
    const suffix = record.status === "template" ? " — disabled: properties unresolved" : "";
    return `<option value="${record.material_id}"${record.material_id===selected?" selected":""}${disabled}>${record.label}${suffix}</option>`;
  }).join("");
}

function material(token) {
  return materialRecords().find((record) => record.material_id===token ||
    (record.aliases || []).includes(token));
}

function esc(text) {
  return String(text).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll('"',"&quot;");
}

function bindInput(element, setter, rerender=false) {
  element.addEventListener("change", () => {
    const value = element.type === "number" ? Number(element.value) : element.value;
    setter(value);
    if (rerender) renderEditors();
    syncJson();
  });
}

function layerDerived(layer) {
  const record = material(layer.material);
  if (!record || !record.density) return {density:"—", areal:"—", geometric:"—", uncertainty:"not available"};
  const density = record.density.value;
  let areal, geometric;
  if (layer.thickness.type === "geometric") {
    geometric = layer.thickness.value * (layer.thickness.unit === "cm" ? 10 : 1);
    areal = density * geometric / 10;
  } else {
    areal = layer.thickness.value;
    geometric = 10 * areal / density;
  }
  const u = record.density.uncertainty || {status:"not_provided"};
  const uncertainty = u.status === "not_provided" ? "not provided" : `${u.value} ${u.unit} (${u.type})`;
  return {density:density.toPrecision(5), areal:areal.toPrecision(5),
          geometric:geometric.toPrecision(5), uncertainty};
}

function renderProfiles() {
  $("profile").innerHTML = capabilities.profiles.map((profile) =>
    `<option value="${profile.profile_id}"${profile.profile_id===model.profile_id?" selected":""}>${profile.label}</option>`).join("");
  const profile = capabilities.profiles.find((item) => item.profile_id===model.profile_id);
  $("quantities").innerHTML = "<strong>Computed quantities:</strong> " + profile.allowed_quantities.map((quantity) =>
    `<label><input type="checkbox" data-quantity="${quantity}" ${model.quantities.includes(quantity)?"checked":""}>${quantity}</label>`).join("");
  document.querySelectorAll("[data-quantity]").forEach((input) => input.addEventListener("change", () => {
    model.quantities = [...document.querySelectorAll("[data-quantity]:checked")].map((x) => x.dataset.quantity);
    syncJson();
  }));
  $("profile-reason").textContent = "Disabled combination: quantities from other profiles have different source/response contracts and cannot be merged here.";
}

function renderEditors() {
  $("shields").innerHTML = model.shield_configurations.map((shield, si) => `
    <article class="editor" data-shield="${si}">
      <div class="editor-head"><h3>Shield ${si+1}</h3><button data-remove-shield="${si}">Remove</button></div>
      <div class="grid three">
        <label>Stable ID<input data-s="${si}" data-field="shield_id" value="${esc(shield.shield_id)}"></label>
        <label>Label<input data-s="${si}" data-field="label" value="${esc(shield.label)}"></label>
        <label>Stack size X (mm)<input type="number" data-s="${si}" data-field="size_x" value="${shield.size_x.value}"></label>
        <label>Stack size Y (mm)<input type="number" data-s="${si}" data-field="size_y" value="${shield.size_y.value}"></label>
        <label>Front X (mm)<input type="number" data-s="${si}" data-field="tx" value="${shield.transform.translation.x.value}"></label>
        <label>Front Y (mm)<input type="number" data-s="${si}" data-field="ty" value="${shield.transform.translation.y.value}"></label>
        <label>Front Z (mm)<input type="number" data-s="${si}" data-field="tz" value="${shield.transform.translation.z.value}"></label>
      </div>
      <p class="transform">Rotation quaternion local→world: w=${shield.transform.rotation_quaternion_local_to_world.w}, x=${shield.transform.rotation_quaternion_local_to_world.x}, y=${shield.transform.rotation_quaternion_local_to_world.y}, z=${shield.transform.rotation_quaternion_local_to_world.z}. Edit non-identity rotations in advanced JSON.</p>
      <h3>Layers (source → absorber)</h3>
      ${shield.layers.map((layer, li) => { const d=layerDerived(layer); return `
        <div class="layer">
          <span class="order">${li+1}</span>
          <label>Layer ID<input data-s="${si}" data-l="${li}" data-layer="layer_id" value="${esc(layer.layer_id)}"></label>
          <label>Material<select data-s="${si}" data-l="${li}" data-layer="material">${materialOptions("shield", material(layer.material)?.material_id || layer.material)}</select></label>
          <label>Input type<select data-s="${si}" data-l="${li}" data-layer="type"><option value="geometric" ${layer.thickness.type==="geometric"?"selected":""}>Geometric</option><option value="areal_density" ${layer.thickness.type==="areal_density"?"selected":""}>Mass / areal</option></select></label>
          <label>Input value (${layer.thickness.type==="geometric"?"mm":"g/cm²"})<input type="number" min="0" step="any" data-s="${si}" data-l="${li}" data-layer="value" value="${layer.thickness.value}"></label>
          <button data-remove-layer="${si},${li}">Remove</button>
          <span></span><small>ρ ${d.density} g/cm³<br>u(ρ): ${d.uncertainty}</small><small>Mass thickness<br>${d.areal} g/cm²</small><small>Geometric<br>${d.geometric} mm</small>
        </div>`; }).join("")}
      <button data-add-layer="${si}" ${shield.layers.length>=5?"disabled title='Maximum five layers'":""}>Add layer</button>
    </article>`).join("");

  $("absorbers").innerHTML = model.absorbers.map((absorber, ai) => `
    <article class="editor" data-absorber="${ai}">
      <div class="editor-head"><h3>Absorber ${ai+1} — slab provider</h3><button data-remove-absorber="${ai}">Remove</button></div>
      <div class="grid absorber-grid">
        <label>Stable ID<input data-a="${ai}" data-field="absorber_id" value="${esc(absorber.absorber_id)}"></label>
        <label>Label<input data-a="${ai}" data-field="label" value="${esc(absorber.label)}"></label>
        <label>Material<select data-a="${ai}" data-field="material">${materialOptions("absorber", material(absorber.material)?.material_id || absorber.material)}</select></label>
        <label>Size X (mm)<input type="number" min="0" data-a="${ai}" data-field="size_x" value="${absorber.geometry.size_x.value}"></label>
        <label>Size Y (mm)<input type="number" min="0" data-a="${ai}" data-field="size_y" value="${absorber.geometry.size_y.value}"></label>
        <label>Thickness (mm)<input type="number" min="0" data-a="${ai}" data-field="thickness" value="${absorber.geometry.thickness.value}"></label>
        <label>Front X (mm)<input type="number" data-a="${ai}" data-field="tx" value="${absorber.geometry.transform.translation.x.value}"></label>
        <label>Front Y (mm)<input type="number" data-a="${ai}" data-field="ty" value="${absorber.geometry.transform.translation.y.value}"></label>
        <label>Front Z (mm)<input type="number" data-a="${ai}" data-field="tz" value="${absorber.geometry.transform.translation.z.value}"></label>
      </div>
      <p class="transform">Provider: slab. Rotation quaternion local→world: identity in the form; advanced JSON accepts an explicit unit quaternion.</p>
    </article>`).join("");
  bindEditorEvents();
}

function bindEditorEvents() {
  document.querySelectorAll("[data-s][data-field]").forEach((input) => bindInput(input, (value) => {
    const shield=model.shield_configurations[Number(input.dataset.s)], field=input.dataset.field;
    if (["size_x","size_y"].includes(field)) shield[field]={value,unit:"mm"};
    else if (["tx","ty","tz"].includes(field)) shield.transform.translation[field[1]]={value,unit:"mm"};
    else shield[field]=value;
  }));
  document.querySelectorAll("[data-layer]").forEach((input) => bindInput(input, (value) => {
    const layer=model.shield_configurations[+input.dataset.s].layers[+input.dataset.l];
    if (input.dataset.layer==="type") { layer.thickness.type=value; layer.thickness.unit=value==="geometric"?"mm":"g/cm2"; }
    else if (input.dataset.layer==="value") layer.thickness.value=value;
    else layer[input.dataset.layer]=value;
  }, true));
  document.querySelectorAll("[data-a][data-field]").forEach((input) => bindInput(input, (value) => {
    const absorber=model.absorbers[+input.dataset.a], field=input.dataset.field;
    if (["size_x","size_y","thickness"].includes(field)) absorber.geometry[field]={value,unit:"mm"};
    else if (["tx","ty","tz"].includes(field)) absorber.geometry.transform.translation[field[1]]={value,unit:"mm"};
    else absorber[field]=value;
  }));
  document.querySelectorAll("[data-remove-shield]").forEach((button) => button.onclick=()=>{ model.shield_configurations.splice(+button.dataset.removeShield,1); renderEditors(); syncJson(); });
  document.querySelectorAll("[data-remove-absorber]").forEach((button) => button.onclick=()=>{ model.absorbers.splice(+button.dataset.removeAbsorber,1); renderEditors(); syncJson(); });
  document.querySelectorAll("[data-remove-layer]").forEach((button) => button.onclick=()=>{ const [s,l]=button.dataset.removeLayer.split(",").map(Number); model.shield_configurations[s].layers.splice(l,1); renderEditors(); syncJson(); });
  document.querySelectorAll("[data-add-layer]").forEach((button) => button.onclick=()=>{ const s=+button.dataset.addLayer; model.shield_configurations[s].layers.push({layer_id:`layer_${model.shield_configurations[s].layers.length+1}`,material:"Al_natural_G4_11_4",thickness:{type:"geometric",value:1,unit:"mm"}}); renderEditors(); syncJson(); });
}

function syncTop() {
  model.request_id=$("request-id").value; model.purpose=$("purpose").value;
  model.accuracy_target=Number($("accuracy").value);
}
function syncJson() { syncTop(); $("request-json").value=JSON.stringify(model,null,2); }
function renderAll() {
  $("request-id").value=model.request_id; $("purpose").value=model.purpose;
  $("accuracy").value=model.accuracy_target ?? 0.10;
  renderProfiles(); renderEditors(); syncJson();
}
function summary(value) {
  if (value.case_matrix) {
    const m=value.case_matrix, c=value.cost_preflight;
    const badges=value.validation_badges.map((b)=>`${b.quantity}: ${b.status} — ${b.reason}`).join("\n");
    return `${m.shield_count} shields × ${m.absorber_count} absorbers = ${m.transport_case_count} transport cases\n`+
      `${c.requested_quantity_series_count} requested quantity series; ${c.minimum_primary_histories}–${c.maximum_primary_histories} histories\n`+
      `Wall time: ${c.wall_time_status}\nManifest: ${value.manifest_content_hash}\n\nValidation badges:\n${badges}\n\n`+
      JSON.stringify(value.uncertainty_preflight,null,2);
  }
  return JSON.stringify(value,null,2);
}

async function act(path, extra={}) {
  try { syncTop(); const value=await api(path,{request:model,bundle_name:$("bundle-name").value,...extra}); $("summary").textContent=summary(value); return value; }
  catch(error) { $("summary").textContent=`ERROR: ${error.message}`; }
}

async function initialize() {
  [capabilities,materialLibrary,model,serverCapability]=await Promise.all([
    api("/api/capabilities"), api("/api/materials"), api("/api/example"),
    api("/api/server")]);
  renderAll();
  if (!serverCapability.submission_enabled) {
    for (const id of ["run","resume","cancel","status"]) {
      $(id).disabled=true; $(id).title=serverCapability.submission_disabled_reason;
    }
    $("profile-reason").textContent += ` Submission disabled: ${serverCapability.submission_disabled_reason}`;
  }
  $("profile").onchange=()=>{ model.profile_id=$("profile").value; model.quantities=[...capabilities.profiles.find((p)=>p.profile_id===model.profile_id).allowed_quantities]; renderProfiles(); syncJson(); };
  ["request-id","purpose","accuracy"].forEach((id)=>$(id).onchange=syncJson);
  $("add-shield").onclick=()=>{ model.shield_configurations.push({shield_id:`shield_${model.shield_configurations.length+1}`,label:"New shield",size_x:mm(50),size_y:mm(50),transform:identity(-20),layers:[{layer_id:"layer_1",material:"Al_natural_G4_11_4",thickness:{type:"geometric",value:1,unit:"mm"}}]}); renderEditors(); syncJson(); };
  $("add-absorber").onclick=()=>{ model.absorbers.push({absorber_id:`absorber_${model.absorbers.length+1}`,label:"New slab",material:"Si_natural_G4_11_4",geometry:{type:"slab",size_x:mm(50),size_y:mm(50),thickness:mm(1),transform:identity(-16)}}); renderEditors(); syncJson(); };
  $("preview").onclick=()=>act("/api/preview");
  $("compile").onclick=()=>act("/api/compile");
  $("run").onclick=()=>act("/api/job/start",{resume:false});
  $("cancel").onclick=()=>act("/api/job/cancel");
  $("resume").onclick=()=>act("/api/job/start",{resume:true});
  $("status").onclick=async()=>{ try { const v=await api(`/api/job?bundle_name=${encodeURIComponent($("bundle-name").value)}`); $("summary").textContent=summary(v); } catch(e) { $("summary").textContent=`ERROR: ${e.message}`; } };
  $("load-json").onclick=()=>{ try { model=JSON.parse($("request-json").value); renderAll(); } catch(e) { $("summary").textContent=`ERROR: invalid request JSON: ${e.message}`; } };
}

initialize().catch((error)=>{ $("summary").textContent=`Initialization failed: ${error.message}`; });
