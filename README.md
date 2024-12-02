Pipeline to check for an enrichment of variants with AF < 0.01 and CADD / splicAI scores > 20 / 0.2 in case cohort in comparison to gnomAD reference cohort.
This is beneficial if the case cohort is small and therefore statistical power is low to detect variants in disease causing genes. Using this approach, caucasian genome samples from gnomAD (~68000) are integrated as a control cohort.

Individual sample level information is lost though. This leads to the downside that multiple variants in the same individual might artificially inflate differences. Linkage disquilibrium can also not be accounted for.
