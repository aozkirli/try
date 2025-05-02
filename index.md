
# Rethinking serial dependence: A large-scale analysis of its effects on the variability of perceptual estimates  


Welcome to the **Rethinking serial dependence: A large-scale analysis of its effects on the variability of perceptual estimates** repository. This project allows users to analyze the effects of serial dependence on the response variability in adjustment tasks

## Main Results (PDF)

[View the main results PDF](figures/main_results_polyfit.pdf)

## Summary of Studies

The table below is generated from `tables/summary_studies.csv`:

<table id="summary-table">
  <thead></thead>
  <tbody></tbody>
</table>
<script>
  document.addEventListener('DOMContentLoaded', () => {
    fetch('tables/summary_studies.csv')
      .then(response => {
        if (!response.ok) throw new Error('Network response was not ok');
        return response.text();
      })
      .then(text => {
        const table = document.getElementById('summary-table');
        const lines = text.trim().split('
');
        if (lines.length === 0) return;
        // Parse header
        const headers = lines[0].split(',');
        const thead = table.querySelector('thead');
        const headerRow = document.createElement('tr');
        headers.forEach(h => {
          const th = document.createElement('th'); th.textContent = h; headerRow.appendChild(th);
        });
        thead.appendChild(headerRow);
        // Parse body rows
        const tbody = table.querySelector('tbody');
        lines.slice(1).forEach(line => {
          const row = document.createElement('tr');
          line.split(',').forEach(cell => {
            const td = document.createElement('td'); td.textContent = cell; row.appendChild(td);
          });
          tbody.appendChild(row);
        });
      })
      .catch(err => console.error('Error loading CSV:', err));
  });
</script>

## Figures Gallery

![Figure 1: Response scatter as a function of feature distance](figures/main_results_polyfit.png)
