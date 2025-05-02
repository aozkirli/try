---
title: "Project Title"
layout: default
---

# Project Title

Welcome to the **Project Title** repository. This project processes user-uploaded datasets through a suite of MATLAB functions and automatically updates result figures and summaries.

## Main Results (PDF)

[View the main results PDF](figures/main_results_polyfit.pdf)

## Summary of Studies

The table below is generated from `tables/summary_studies.csv`:

```html
<div id="summary-table"></div>
<script>
  fetch('tables/summary_studies.csv')
    .then(response => response.text())
    .then(text => {
      const table = document.getElementById('summary-table');
      const rows = text.trim().split('\n');
      rows.forEach((row, rowIndex) => {
        const tr = document.createElement('tr');
        row.split(',').forEach(cell => {
          const cellElem = document.createElement(rowIndex === 0 ? 'th' : 'td');
          cellElem.textContent = cell;
          tr.appendChild(cellElem);
        });
        table.appendChild(tr);
      });
    })
    .catch(err => console.error('Error loading CSV:', err));
</script>
```

## Figures Gallery

![Figure 1: Example Description](figures/figure1.png)

## Automating Updates

See `.github/workflows/update-index.yml` and `scripts/generate_index.py` for details on automatically regenerating this file when figures or tables change.
