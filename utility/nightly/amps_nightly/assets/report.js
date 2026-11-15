/* Status filter buttons on family pages; the pages work without this script. */
(function () {
  'use strict';
  function init() {
    var bar = document.querySelector('[data-filters]');
    if (!bar) return;
    var buttons = Array.prototype.slice.call(bar.querySelectorAll('[data-filter]'));
    var rows = Array.prototype.slice.call(document.querySelectorAll('[data-status]'));
    var empty = document.querySelector('[data-empty]');
    function apply(filter) {
      var shown = 0;
      rows.forEach(function (row) {
        var visible = filter === 'all' || row.getAttribute('data-status') === filter;
        row.hidden = !visible;
        if (visible) shown++;
      });
      buttons.forEach(function (b) { b.setAttribute('aria-pressed', String(b.getAttribute('data-filter') === filter)); });
      if (empty) empty.hidden = shown > 0;
      try { history.replaceState(null, '', filter === 'all' ? location.pathname : '#' + filter); } catch (e) { /* file:// */ }
    }
    buttons.forEach(function (b) {
      b.addEventListener('click', function () { apply(b.getAttribute('data-filter')); });
    });
    bar.hidden = false;
    var initial = location.hash.replace('#', '');
    apply(buttons.some(function (b) { return b.getAttribute('data-filter') === initial; }) ? initial : 'all');
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init); else init();
})();
