# =====================================================================
# Photographic specimen catalogue  —  Sparidae field season 2026
# =====================================================================
# Builds ONE self-contained, searchable/filterable HTML file from the same
# data analysis.R uses, plus the photos in all_records/photos/.
#
# RE-RUN WORKFLOW
#   Add new bream_*.csv exports to all_records/ and the matching photos to
#   all_records/photos/, then source this file again. It re-globs everything
#   and rebuilds catalog/fish_catalog.html from scratch.
#
# PHOTO MATCHING
#   Photos are located by SCANNING all_records/photos/ and matching each file
#   back to a fish_id, after stripping a trailing " copy" / " copy 2" and the
#   extension. So duplicated files (e.g. "<id> copy.jpg") and files with the
#   extension hidden/removed all still resolve.
#
# REQUIRES:  install.packages(c("magick", "jsonlite"))   # plus tidyverse/here/janitor
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(janitor)
})

if (!requireNamespace("magick", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("This script needs the 'magick' and 'jsonlite' packages:\n",
       "  install.packages(c(\"magick\", \"jsonlite\"))", call. = FALSE)
}

# ---- config (tune to trade quality vs. file size) -------------------
MAX_PX   <- 760L   # longest edge of each embedded photo (lower -> smaller HTML)
QUALITY  <- 80L    # embedded JPEG quality
OUT_HTML <- here("catalog", "fish_catalog.html")

# Cross-dataset de-duplication, matching 02_analysis.R. A fish caught during the season and also
# logged by the government tagging programme carries a different id in each file, so it survives
# distinct(fish_id) as two records. 01_combine_data.R resolves those and 02_analysis.R now honours
# that resolution; the catalogue must too, or its specimen count contradicts the paper's.
#
# THE TRADE-OFF, stated because it is not free: the resolution keeps the historical row, which has
# no photograph, and discards the app row, which does. Filtering therefore costs the catalogue one
# photograph. Set FALSE to keep every photographed record, and expect the catalogue's n to exceed
# the paper's by the number of duplicates.
DEDUPE_CATALOG <- TRUE
COMBINED_CSV   <- here("combined_records", "combined_dataset.csv")
SEASON_LABEL   <- "2026 season"

# Vernacular names: the HMGoG list is canonical and the app uses angler-facing names. Same map as
# 01_combine_data.R and 02_analysis.R, so the catalogue, the figures and the tables cannot label
# one species three ways. Keep these three copies identical.
common_fix <- c(
  "Pagrus auriga"              = "Redbanded Seabream",
  "Diplodus cervinus cervinus" = "Soldier Seabream",
  "Pagellus acarne"            = "Bronze Seabream",
  "Pagrus pagrus"              = "Common Seabream"
)

# ---- 1. LOAD  (same auto-discovering loader as analysis.R) ----------
raw <- list.files(here("all_records"), pattern = "^bream_.*\\.csv$", full.names = TRUE) |>
  set_names(basename) |>
  map(\(f) read_csv(f, na = c("", "NA", "NR"), show_col_types = FALSE)) |>
  list_rbind() |>
  clean_names() |>
  distinct(fish_id, .keep_all = TRUE) |>
  mutate(common_name = coalesce(unname(common_fix[scientific_name]), common_name))

if (isTRUE(DEDUPE_CATALOG)) {
  if (file.exists(COMBINED_CSV)) {
    keep <- read_csv(COMBINED_CSV, na = c("", "NA", "NR"), show_col_types = FALSE) |>
      clean_names() |>
      dplyr::filter(dataset == SEASON_LABEL, !is.na(fish_id)) |>
      dplyr::pull(fish_id)
    dropped <- setdiff(raw$fish_id, keep)
    raw <- raw |> dplyr::filter(fish_id %in% keep)
    if (length(dropped))
      message(sprintf("[catalog] %d cross-dataset duplicate(s) removed: %s. Catalogue n = %d, matching the analysis.",
                      length(dropped), paste(dropped, collapse = ", "), nrow(raw)))
  } else {
    warning("[catalog] ", basename(COMBINED_CSV), " not found - de-duplication SKIPPED. ",
            "The catalogue's specimen count may exceed the analysis's. Run 01_combine_data.R first.",
            call. = FALSE)
  }
}

# ---- 2. ASSEMBLE SPECIMEN RECORDS -----------------------------------
catalog <- raw |>
  mutate(
    date_d = ymd(date),
    length_label = case_when(
      !is.na(length_true_cm)  ~ paste0(length_true_cm, " cm"),
      !is.na(length_estimate) ~ paste0("~", length_estimate, " cm (est.)"),
      TRUE                    ~ NA_character_
    ),
    # numeric length for sorting: measured value, or midpoint of an estimate range
    length_sort = coalesce(
      length_true_cm,
      map_dbl(str_split(length_estimate, "-"), \(x) suppressWarnings(mean(as.numeric(x))))
    ),
    length_sort  = ifelse(is.nan(length_sort), NA_real_, length_sort),
    weight_label = if_else(!is.na(weight_kg), paste0(weight_kg, " kg"), NA_character_),
    temp_label   = if_else(!is.na(water_temp_c), paste0(water_temp_c, " \u00b0C"), NA_character_),
    loc          = if_else(!is.na(latitude) & !is.na(longitude),
                           sprintf("%.4f, %.4f", latitude, longitude), NA_character_),
    release_fate = as.character(release_fate),
    tag_status = case_when(
      as.logical(was_tagged)       ~ "Newly tagged",
      as.logical(had_existing_tag) ~ "Recapture",
      TRUE                         ~ "Untagged"
    ),
    tag_id     = coalesce(new_tag_id, existing_tag_id),
    date_iso   = as.character(date_d),
    date_label = format(date_d, "%d %b %Y")
  )

# ---- 3. MATCH PHOTOS ON DISK  (robust to " copy" / missing extension) ----
photo_dir <- here("all_records", "photos")
photo_index <- tibble(path = list.files(photo_dir, full.names = TRUE)) |>
  filter(!dir.exists(path), !str_detect(basename(path), "^\\.")) |>  # skip subfolders & hidden files
  mutate(
    key = path |>
      basename() |>
      tools::file_path_sans_ext() |>
      str_remove(regex("\\s*copy\\s*\\d*\\s*$", ignore_case = TRUE)) |>  # drop a " copy N" suffix
      str_remove(regex("_\\d+$")) |>                                     # drop a _1 / _2 multi-photo suffix
      str_squish()
  ) |>
  filter(key != "") |>
  arrange(path) |>
  group_by(key) |>
  summarise(photo_paths = list(path), .groups = "drop")   # ALL photos per fish, in filename order

catalog <- catalog |>
  left_join(photo_index |> rename(fish_id = key), by = "fish_id")

# ---- 4. EMBED PHOTOS as resized data-URI thumbnails -----------------
encode_photo <- function(p) {
  if (is.na(p) || !file.exists(p)) return(NA_character_)
  tryCatch({
    img <- magick::image_read(p)
    img <- magick::image_scale(img, paste0(MAX_PX, "x", MAX_PX, ">"))  # shrink only
    buf <- magick::image_write(img, format = "jpeg", quality = QUALITY)
    paste0("data:image/jpeg;base64,", jsonlite::base64_enc(buf))
  }, error = function(e) NA_character_)
}

message("Matched ", sum(lengths(catalog$photo_paths) > 0), " specimens with photos; encoding\u2026")
catalog <- catalog |>
  mutate(imgs = map(photo_paths, function(ps) {
    if (is.null(ps) || !length(ps)) return(character(0))
    enc <- vapply(ps, encode_photo, character(1))
    unname(enc[!is.na(enc)])                       # keep every image that encoded cleanly
  }))
n_embedded <- sum(lengths(catalog$imgs) > 0)
message("Embedded photos for ", n_embedded, " of ", nrow(catalog), " specimens (",
        sum(lengths(catalog$imgs)), " images total).")
missing_ids <- catalog |> filter(lengths(imgs) == 0) |> pull(fish_id)
if (length(missing_ids))
  message("No photo for: ", paste(head(missing_ids, 20), collapse = ", "),
          if (length(missing_ids) > 20) " \u2026" else "")

# ---- 5. RECORDS + STATS as JSON -------------------------------------
records <- catalog |>
  arrange(desc(date_iso), common_name) |>
  transmute(
    fish_id, common_name, scientific_name, spanish_name, genus,
    length_label, length_sort, weight_label,
    date_iso, date_label, time, loc, temp_label,
    release_fate, tag_status, tag_id, imgs
  )
data_json <- jsonlite::toJSON(records, na = "null", auto_unbox = TRUE)

stats <- list(
  n            = nrow(catalog),
  species      = n_distinct(catalog$common_name),
  photographed = n_embedded,
  date_min     = format(min(catalog$date_d, na.rm = TRUE), "%d %b %Y"),
  date_max     = format(max(catalog$date_d, na.rm = TRUE), "%d %b %Y")
)
stats_json <- jsonlite::toJSON(stats, auto_unbox = TRUE)

# ---- 6. ASSEMBLE + WRITE SELF-CONTAINED HTML ------------------------
PAGE_TOP <- r"---(<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>The Sparidae Catalogue &mdash; Field Season 2026</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght@0,9..144,400;0,9..144,600;0,9..144,700;1,9..144,400;1,9..144,500&family=Hanken+Grotesk:wght@400;500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
:root{
  --paper:#efe7d8; --paper-2:#e6dcc7; --paper-3:#dccfb6;
  --card:#fbf7ef; --ink:#17120b; --ink-soft:#5b5345; --ink-faint:#8d836f;
  --teal:#1f7077; --teal-deep:#124a4f; --ochre:#b8923f;
  --line:rgba(23,18,11,.13); --line-2:rgba(23,18,11,.22);
  --shadow:0 1px 2px rgba(23,18,11,.06), 0 8px 22px rgba(23,18,11,.05);
  --shadow-hi:0 5px 12px rgba(23,18,11,.10), 0 20px 52px rgba(23,18,11,.13);
  --font-display:"Fraunces","Georgia",serif;
  --font-body:"Hanken Grotesk",system-ui,sans-serif;
  --font-mono:"IBM Plex Mono",ui-monospace,monospace;
}
*{box-sizing:border-box}
html{scroll-behavior:smooth}
body{margin:0;color:var(--ink);font-family:var(--font-body);line-height:1.45;
  background:radial-gradient(135% 90% at 50% -12%,var(--paper) 0%,var(--paper-2) 64%,var(--paper-3) 100%) fixed;
  -webkit-font-smoothing:antialiased;text-rendering:optimizeLegibility;}
body::before{content:"";position:fixed;inset:0;z-index:0;pointer-events:none;opacity:.05;mix-blend-mode:multiply;
  background-image:url("data:image/svg+xml,%3Csvg%20xmlns='http://www.w3.org/2000/svg'%20width='180'%20height='180'%3E%3Cfilter%20id='n'%3E%3CfeTurbulence%20type='fractalNoise'%20baseFrequency='0.8'%20numOctaves='2'%20stitchTiles='stitch'/%3E%3C/filter%3E%3Crect%20width='100%25'%20height='100%25'%20filter='url(%23n)'/%3E%3C/svg%3E");}
.wrap{position:relative;z-index:1;max-width:1260px;margin:0 auto;padding:0 26px 90px;}

header.masthead{padding:58px 0 26px;}
.kicker{font-family:var(--font-mono);font-size:.72rem;letter-spacing:.3em;text-transform:uppercase;color:var(--teal-deep);}
.masthead h1{font-family:var(--font-display);font-weight:600;font-size:clamp(2.3rem,5.4vw,3.9rem);line-height:1.0;margin:.2em 0 .12em;letter-spacing:-.015em;}
.masthead h1 em{font-style:italic;font-weight:500;color:var(--teal-deep);}
.masthead .sub{max-width:62ch;color:var(--ink-soft);font-size:1.02rem;}
.stats{display:flex;flex-wrap:wrap;gap:34px;margin-top:26px;padding-top:22px;border-top:1px solid var(--line-2);}
.stat .num{font-family:var(--font-display);font-weight:600;font-size:1.75rem;line-height:1;}
.stat.wide .num{font-family:var(--font-mono);font-weight:500;font-size:1rem;letter-spacing:.01em;padding-top:.5em;}
.stat .lab{font-family:var(--font-mono);font-size:.64rem;letter-spacing:.2em;text-transform:uppercase;color:var(--ink-faint);margin-top:7px;}

.toolbar{position:sticky;top:0;z-index:30;display:flex;flex-wrap:wrap;gap:11px;align-items:center;
  padding:14px 0;margin-bottom:28px;border-bottom:1px solid var(--line);
  background:linear-gradient(var(--paper),rgba(239,231,216,.85));backdrop-filter:blur(9px);}
.toolbar .search{flex:1 1 230px;}
.toolbar input[type=search]{width:100%;font-family:var(--font-body);font-size:.95rem;padding:11px 14px;color:var(--ink);
  background:var(--card);border:1px solid var(--line-2);border-radius:3px;}
.toolbar input[type=search]::placeholder{color:var(--ink-faint);}
.toolbar select{font-family:var(--font-body);font-size:.86rem;padding:10px 12px;color:var(--ink);background:var(--card);
  border:1px solid var(--line-2);border-radius:3px;cursor:pointer;}
.chk{font-family:var(--font-mono);font-size:.68rem;letter-spacing:.08em;text-transform:uppercase;color:var(--ink-soft);
  display:flex;align-items:center;gap:7px;cursor:pointer;user-select:none;}
.chk input{accent-color:var(--teal);}
.count{margin-left:auto;font-family:var(--font-mono);font-size:.72rem;letter-spacing:.14em;text-transform:uppercase;color:var(--ink-faint);}
input:focus,select:focus{outline:2px solid var(--teal);outline-offset:1px;}

.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:22px;}
.card{position:relative;background:var(--card);border:1px solid var(--line);border-radius:5px;overflow:hidden;cursor:pointer;
  box-shadow:var(--shadow);transition:transform .26s cubic-bezier(.2,.7,.2,1),box-shadow .26s,border-color .26s;
  animation:rise .5s ease both;animation-delay:calc(var(--i,0)*32ms);}
.card:hover,.card:focus-visible{transform:translateY(-4px);box-shadow:var(--shadow-hi);border-color:var(--teal);outline:none;}
.ph{aspect-ratio:4/3;background:#e7ddca;overflow:hidden;position:relative;}
.ph img{width:100%;height:100%;object-fit:cover;display:block;transition:transform .45s ease;}
.card:hover .ph img{transform:scale(1.045);}
.ph-empty{display:flex;align-items:center;justify-content:center;color:var(--ink-faint);
  background:repeating-linear-gradient(45deg,#e7ddca,#e7ddca 11px,#e2d7c1 11px,#e2d7c1 22px);}
.ph-empty span{font-family:var(--font-mono);font-size:.66rem;letter-spacing:.16em;text-transform:uppercase;}
.ph-count{position:absolute;bottom:8px;right:8px;background:rgba(30,30,30,.72);color:#fff;font-family:var(--font-mono);font-size:.64rem;letter-spacing:.06em;padding:2px 7px;border-radius:10px;}
.body{padding:14px 15px 12px;}
.cn{font-family:var(--font-display);font-weight:600;font-size:1.13rem;line-height:1.12;margin:0;}
.sn{font-family:var(--font-display);font-style:italic;font-weight:400;color:var(--ink-soft);font-size:.92rem;margin:.16em 0 .72em;}
dl.data{margin:0 0 11px;display:grid;gap:3px;}
.data .row{display:flex;gap:9px;font-size:.82rem;}
.data dt{font-family:var(--font-mono);font-size:.6rem;letter-spacing:.13em;text-transform:uppercase;color:var(--ink-faint);width:46px;flex:none;padding-top:3px;}
.data dd{margin:0;color:var(--ink);}
.badges{display:flex;flex-wrap:wrap;gap:6px;margin-bottom:11px;}
.badge,.tag{font-family:var(--font-mono);font-size:.6rem;letter-spacing:.1em;text-transform:uppercase;padding:3px 8px;border-radius:2px;}
.badge.kept{background:var(--ochre);color:#2a1f06;}
.badge.released{background:transparent;color:var(--teal-deep);border:1px solid var(--teal);}
.tag{background:rgba(23,18,11,.05);color:var(--ink-soft);}
.tag.newtag{color:var(--teal-deep);background:rgba(31,112,119,.1);}
.tag.recap{color:#8a5a12;background:rgba(184,146,63,.16);}
.acc{font-family:var(--font-mono);font-size:.63rem;letter-spacing:.06em;color:var(--ink-faint);padding-top:9px;border-top:1px solid var(--line);word-break:break-all;}
.empty{grid-column:1/-1;text-align:center;padding:74px 0;font-family:var(--font-display);font-style:italic;font-size:1.35rem;color:var(--ink-soft);}

.lb{position:fixed;inset:0;z-index:100;display:none;align-items:center;justify-content:center;padding:30px;
  background:rgba(20,15,9,.74);backdrop-filter:blur(3px);}
.lb.open{display:flex;animation:fade .25s ease;}
.lb-panel{position:relative;max-width:960px;width:100%;max-height:88vh;border-radius:6px;overflow:hidden;
  box-shadow:0 30px 84px rgba(0,0,0,.42);animation:pop .3s cubic-bezier(.2,.8,.2,1);}
.lb-body{display:grid;grid-template-columns:1.15fr 1fr;max-height:88vh;overflow:auto;background:var(--card);}
.lb-img{background:#1b1610;display:flex;align-items:center;justify-content:center;}
.lb-img img{width:100%;height:100%;max-height:88vh;object-fit:contain;display:block;}
.lb-meta{padding:32px 32px 30px;}
.lb-meta h2{font-family:var(--font-display);font-weight:600;font-size:1.8rem;margin:0;line-height:1.05;}
.lb-meta .sn{font-size:1.06rem;margin:.22em 0 1.15em;}
.lb-meta dl{display:grid;margin:0 0 18px;}
.lrow{display:flex;gap:14px;padding:8px 0;border-bottom:1px solid var(--line);}
.lrow dt{font-family:var(--font-mono);font-size:.62rem;letter-spacing:.13em;text-transform:uppercase;color:var(--ink-faint);width:112px;flex:none;padding-top:4px;}
.lrow dd{margin:0;}
.lb .acc{border:none;padding:0 0 16px;}
.dl{display:inline-block;font-family:var(--font-mono);font-size:.7rem;letter-spacing:.12em;text-transform:uppercase;text-decoration:none;
  color:var(--paper);background:var(--teal-deep);padding:10px 17px;border-radius:3px;transition:background .2s;}
.dl:hover{background:var(--teal);}
.lb-close{position:absolute;top:12px;right:14px;z-index:3;width:34px;height:34px;border:none;border-radius:50%;cursor:pointer;
  background:rgba(255,255,255,.92);font-size:21px;line-height:1;color:#222;}
.lb-close:hover{background:#fff;}
@media (max-width:720px){.lb-body{grid-template-columns:1fr;}.lb-img img{max-height:46vh;}}

@keyframes rise{from{opacity:0;transform:translateY(12px);}to{opacity:1;transform:translateY(0);}}
@keyframes fade{from{opacity:0;}to{opacity:1;}}
@keyframes pop{from{opacity:0;transform:scale(.97) translateY(8px);}to{opacity:1;transform:none;}}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important;}}
</style>
</head>
<body>
<div class="wrap">
  <header class="masthead">
    <div class="kicker">Field Season 2026 &middot; Bay of Gibraltar</div>
    <h1>The Sparidae <em>Catalogue</em></h1>
    <p class="sub">A working photographic record of the sea bream sampled during the 2026 field season. Each specimen is keyed to its unique capture ID, with morphometric, locational and tagging data.</p>
    <div class="stats" id="stats"></div>
  </header>
  <div class="toolbar">
    <div class="search"><input id="q" type="search" placeholder="Search species, ID, tag&hellip;" aria-label="Search"></div>
    <select id="f-species" aria-label="Species"><option value="">All species</option></select>
    <select id="f-fate" aria-label="Release fate"><option value="">Any fate</option><option>Kept</option><option>Released Alive</option></select>
    <select id="f-tag" aria-label="Tag status"><option value="">Any tag status</option><option>Newly tagged</option><option>Recapture</option><option>Untagged</option></select>
    <select id="f-sort" aria-label="Sort order"><option value="date_desc">Newest first</option><option value="date_asc">Oldest first</option><option value="len_desc">Largest first</option><option value="len_asc">Smallest first</option><option value="sp_asc">Species A&ndash;Z</option></select>
    <label class="chk"><input id="f-photo" type="checkbox"> Photographed only</label>
    <span class="count" id="count"></span>
  </div>
  <main class="grid" id="grid"></main>
</div>
<div class="lb" id="lb" role="dialog" aria-modal="true" aria-label="Specimen detail">
  <div class="lb-panel">
    <button class="lb-close" aria-label="Close">&times;</button>
    <div class="lb-body" id="lb-body"></div>
  </div>
</div>
<script>
)---"

PAGE_BOTTOM <- r"---(
const $ = id => document.getElementById(id);
const grid=$("grid"), q=$("q"), fSpecies=$("f-species"), fFate=$("f-fate"),
      fTag=$("f-tag"), fSort=$("f-sort"), fPhoto=$("f-photo"), countEl=$("count");

function esc(s){return (s==null?"":String(s)).replace(/[&<>"]/g,c=>({"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;"}[c]));}

// header stats
const statItems=[[STATS.n,"Specimens"],[STATS.species,"Species"],[STATS.photographed,"Photographed"],
  [STATS.date_min+" \u2013 "+STATS.date_max,"Sampling window"]];
$("stats").innerHTML=statItems.map((s,i)=>`<div class="stat${i===3?" wide":""}"><div class="num">${esc(s[0])}</div><div class="lab">${s[1]}</div></div>`).join("");

// species filter options
[...new Set(DATA.map(d=>d.common_name))].sort().forEach(s=>{
  const o=document.createElement("option");o.value=s;o.textContent=s;fSpecies.appendChild(o);
});

function fateClass(f){return f==="Kept"?"badge kept":"badge released";}
function lenVal(d){return (d.length_sort==null||isNaN(d.length_sort))?-Infinity:d.length_sort;}
function imgList(d){return Array.isArray(d.imgs)?d.imgs:(d.imgs?[d.imgs]:[]);}

function cardHTML(d,i){
  const _imgs=imgList(d);
  const photo=_imgs.length
    ? `<div class="ph"><img loading="lazy" src="${_imgs[0]}" alt="${esc(d.common_name)} \u2014 ${esc(d.fish_id)}">${_imgs.length>1?`<span class="ph-count">${_imgs.length}\u00d7</span>`:``}</div>`
    : `<div class="ph ph-empty"><span>No photograph</span></div>`;
  const tagLine=d.tag_status==="Untagged"
    ? `<span class="tag">Untagged</span>`
    : `<span class="tag ${d.tag_status==="Recapture"?"recap":"newtag"}">${esc(d.tag_status)}${d.tag_id?" \u00b7 "+esc(d.tag_id):""}</span>`;
  const rows=[["Length",d.length_label],["Date",d.date_label],["Temp",d.temp_label]]
    .filter(r=>r[1]).map(r=>`<div class="row"><dt>${r[0]}</dt><dd>${esc(r[1])}</dd></div>`).join("");
  return `<article class="card" tabindex="0" data-id="${esc(d.fish_id)}" style="--i:${Math.min(i,28)}">
    ${photo}
    <div class="body">
      <h3 class="cn">${esc(d.common_name)}</h3>
      <p class="sn">${esc(d.scientific_name)}</p>
      <dl class="data">${rows}</dl>
      <div class="badges"><span class="${fateClass(d.release_fate)}">${esc(d.release_fate)}</span>${tagLine}</div>
      <div class="acc">${esc(d.fish_id)}</div>
    </div></article>`;
}

function current(){
  const term=q.value.trim().toLowerCase(), sp=fSpecies.value, fa=fFate.value, tg=fTag.value, ph=fPhoto.checked;
  let list=DATA.filter(d=>{
    if(sp && d.common_name!==sp) return false;
    if(fa && d.release_fate!==fa) return false;
    if(tg && d.tag_status!==tg) return false;
    if(ph && !imgList(d).length) return false;
    if(term){
      const hay=[d.common_name,d.scientific_name,d.spanish_name,d.genus,d.fish_id,d.tag_id].join(" ").toLowerCase();
      if(!hay.includes(term)) return false;
    }
    return true;
  });
  const s=fSort.value;
  const byDate=(a,b)=>a.date_iso<b.date_iso?1:a.date_iso>b.date_iso?-1:0;
  if(s==="date_desc") list.sort(byDate);
  else if(s==="date_asc") list.sort((a,b)=>-byDate(a,b));
  else if(s==="len_desc") list.sort((a,b)=>lenVal(b)-lenVal(a));
  else if(s==="len_asc") list.sort((a,b)=>lenVal(a)-lenVal(b));
  else if(s==="sp_asc") list.sort((a,b)=>a.common_name.localeCompare(b.common_name));
  return list;
}

function render(){
  const list=current();
  countEl.textContent=list.length+(list.length===1?" specimen":" specimens");
  grid.innerHTML=list.length?list.map(cardHTML).join(""):`<p class="empty">No fish match these filters.</p>`;
}
[q,fSpecies,fFate,fTag,fSort,fPhoto].forEach(el=>el.addEventListener("input",render));

// lightbox
const lb=$("lb"), lbBody=$("lb-body");
function openLB(id){
  const d=DATA.find(x=>x.fish_id===id); if(!d) return;
  const _imgs=imgList(d);
  const img=_imgs.length?_imgs.map(u=>`<img src="${u}" alt="${esc(d.common_name)}" style="display:block;max-width:100%;margin:0 auto 10px;">`).join(""):`<div class="ph-empty" style="min-height:300px"><span>No photograph</span></div>`;
  const rows=[["Scientific name",d.scientific_name],["Spanish name",d.spanish_name],["Genus",d.genus],
    ["Length",d.length_label],["Weight",d.weight_label],["Date",d.date_label],["Time",d.time],
    ["Water temp",d.temp_label],["Position",d.loc],["Release",d.release_fate],
    ["Tag",d.tag_status+(d.tag_id?" \u00b7 "+d.tag_id:"")]]
    .filter(r=>r[1]).map(r=>`<div class="lrow"><dt>${r[0]}</dt><dd>${esc(r[1])}</dd></div>`).join("");
  const dl=imgList(d).length?`<a class="dl" href="${imgList(d)[0]}" download="${esc(d.fish_id)}.jpg">Download image</a>`:"";
  lbBody.innerHTML=`<div class="lb-img">${img}</div>
    <div class="lb-meta"><h2>${esc(d.common_name)}</h2><p class="sn">${esc(d.scientific_name)}</p>
    <dl>${rows}</dl><div class="acc">${esc(d.fish_id)}</div>${dl}</div>`;
  lb.classList.add("open"); document.body.style.overflow="hidden";
}
function closeLB(){lb.classList.remove("open"); document.body.style.overflow="";}
grid.addEventListener("click",e=>{const c=e.target.closest(".card"); if(c) openLB(c.dataset.id);});
grid.addEventListener("keydown",e=>{if((e.key==="Enter"||e.key===" ")&&e.target.classList.contains("card")){e.preventDefault();openLB(e.target.dataset.id);}});
lb.addEventListener("click",e=>{if(e.target===lb||e.target.classList.contains("lb-close")) closeLB();});
document.addEventListener("keydown",e=>{if(e.key==="Escape") closeLB();});

render();
</script>
</body>
</html>
)---"

html <- paste0(PAGE_TOP, "\n",
               "const STATS = ", stats_json, ";\n",
               "const DATA = ",  data_json,  ";\n\n",
               PAGE_BOTTOM)

dir.create(dirname(OUT_HTML), showWarnings = FALSE, recursive = TRUE)
con <- file(OUT_HTML, open = "wb")
writeBin(charToRaw(enc2utf8(html)), con)
close(con)
message("Wrote ", OUT_HTML)
if (interactive()) utils::browseURL(OUT_HTML)