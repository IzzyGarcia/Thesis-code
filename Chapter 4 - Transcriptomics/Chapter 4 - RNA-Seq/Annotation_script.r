---
title: "ZYFS_GFP_vs_Transformed_control_results"
output: html_document
editor_options: 
  chunk_output_type: console
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

#libraries:

```{r}
library(tidyverse)
```

#load in the gene list

```{r}
DE_expression_data <- read.csv("ZYFS_GFP_vs_Transformed_control_results.csv")
```

#read in annotation

```{r}
uniprot_annotation <- read.delim("./uniprot-annotation.tab")

#now we just need to rename the primary gene name column in the uniprot data to match the symbol column in the DE data
#and for convenience, rename the Entry column to Uniprot.ID
uniprot_annotation <- uniprot_annotation %>% rename(SYMBOL = Gene.names...primary..)
uniprot_annotation <- uniprot_annotation %>% rename(Uniprot.ID = Entry)
```

# apply the annotation data to the list (left_join)

```{r}
annotated_DE_expression_data <- left_join(DE_expression_data, uniprot_annotation)
```

# write out

```{r}
write_csv(annotated_DE_expression_data, "./annotated_filtered_ZFYS_GFP_transfected.csv")
