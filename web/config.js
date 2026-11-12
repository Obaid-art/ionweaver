// Where the engine lives.  Served by the engine itself (local run, Hugging Face
// Space): same origin.  Served as a static page (GitHub Pages): the Space's URL.
window.IONWEAVER_API = (function () {
  const h = location.hostname;
  if (h.endsWith("github.io") || location.protocol === "file:") {
    return "https://ionweaver.onrender.com";   // the engine server
  }
  return "";
})();
