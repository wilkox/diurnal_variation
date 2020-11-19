all: seasonality_strengths.tsv compositional_seasonality_strengths.tsv

seasonality_strengths.tsv: calculate_seasonality_strengths.R
	Rscript --vanilla --no-save calculate_seasonality_strengths.R
	
compositional_seasonality_strengths.tsv: calculate_compositional_seasonality_strengths.R
	Rscript --vanilla --no-save calculate_compositional_seasonality_strengths.R

no_nc_seasonality_strengths.tsv: calculate_no_nc_seasonality_strengths.R
	Rscript --vanilla --no-save calculate_no_nc_seasonality_strengths.R
