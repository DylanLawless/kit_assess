library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)

# Import data ----
file="../data/processed/target_info/AH_S1_target_info_biomart_simple.txt"
df = read.csv(file, sep = "\t", )
df <- df %>% arrange(Chromosome.scaffold.name, Gene.start.bp) 	

file="../data/processed/target_info/CH_S2_target_info_biomart_simple.txt"
df2 = read.csv(file, sep = "\t", )
df2 <- df2 %>% arrange(Chromosome.scaffold.name, Gene.start.bp) 	

names(df)

df$AH_S1 <- "Yes"
df2$CH_S2 <- "Yes"

merged <- left_join(df, df2, by=c("Gene.name", "Chromosome.scaffold.name", "Gene.start.bp", "Gene.end.bp", "Gene.pc.GC.content"))


df <- merged
df <- df %>% select(Chromosome.scaffold.name, Gene.start.bp, Gene.end.bp, everything())
df2 <- df2 %>% select(Chromosome.scaffold.name, Gene.start.bp, Gene.end.bp, everything())

#colnames(df)[colnames(df) == 'Database_category'] <- "Database category"
#colnames(df)[colnames(df) == 'Usage'] <- "In use"

#library(stringr)
#df$Usage <- str_replace(df$Usage, "yes", "Yes")




# reactable ----
library(reactable)
options(reactable.theme = reactableTheme(
	borderColor = "#dfe2e5",
	stripedColor = "#E5E5E5",
	highlightColor = "#fcf0e6",
	cellPadding = "8px 12px",
	style = list(fontFamily = "-apple-system, Arial, BlinkMacSystemFont, Segoe UI, Helvetica,  sans-serif",
					 fontSize = "0.8rem"),
	searchInputStyle = list(width = "50%")
))

df_t <- 
	reactable(df,
				  compact = TRUE,
				  searchable = TRUE,
				  elementId = "download-table",
				  defaultPageSize = 170,
				  defaultColDef = colDef(minWidth = 80, maxWidth = 200),
				  columns = list(
				  	"In Use" = colDef(maxWidth = 40) #, 
		#		  	"URL" = colDef(minWidth = 300,
	#			  						cell = function(value, index) {
#				  		url <- sprintf(df[index, "URL"], value)
#				  		htmltools::tags$a(href = url, target = "_blank", as.character(value))
#				  	 # htmltools::tags$a(href = url, target = "_blank", "link")
#				  	}),
#				  	"In use" = colDef(cell = function(value) {
#				  		# Render as an X mark or check mark
#				  		if (value == "No") "\u274c No" else "\u2714\ufe0f Yes"
#				  	})
				  	
				  ),
				 filterable = TRUE,
				 showSortable = TRUE,
				 showPageSizeOptions = TRUE,
				 striped = TRUE,
				 highlight = TRUE
	)

df_t
library(reactablefmtr)
save_reactable(df_t, "../data/processed/target_info/target_info_biomart_simple.html")


write.table(df , file='../data/processed/target_info/target_info_biomart_simple.csv', sep=",", quote=FALSE, row.names=TRUE, col.names = TRUE)


df2_t <- 
	reactable(df2,
				 compact = TRUE,
				 searchable = TRUE,
				 elementId = "download-table",
				 defaultPageSize = 170,
				 defaultColDef = colDef(minWidth = 140),
				 columns = list(
				 	"In Use" = colDef(maxWidth = 40) #, 
				 	#		  	"URL" = colDef(minWidth = 300,
				 	#			  						cell = function(value, index) {
				 	#				  		url <- sprintf(df[index, "URL"], value)
				 	#				  		htmltools::tags$a(href = url, target = "_blank", as.character(value))
				 	#				  	 # htmltools::tags$a(href = url, target = "_blank", "link")
				 	#				  	}),
				 	#				  	"In use" = colDef(cell = function(value) {
				 	#				  		# Render as an X mark or check mark
				 	#				  		if (value == "No") "\u274c No" else "\u2714\ufe0f Yes"
				 	#				  	})
				 	
				 ),
				 filterable = TRUE,
				 showSortable = TRUE,
				 showPageSizeOptions = TRUE,
				 striped = TRUE,
				 highlight = TRUE
	)

df2_t
library(reactablefmtr)
save_reactable(df2_t, "../data/processed/target_info/CH_S2_target_info_biomart_simple.html")



# scratch

file="../data/processed/target_info/AH_S1_target_info_biomart_simple.txt"
df = read.csv(file, sep = "\t", )

df <- unite(df, Chromosome.scaffold.name, Gene.start.bp, col = "Pos", sep = ".")
df$Pos <- as.numeric(df$Pos)
df <- df %>% arrange(Pos) 	

df$value <- rownames(df)
df$value <- as.numeric(df$value)
df <- df %>% select(Pos, value, Gene.name)
df$Region <- "coding"

file="../data/processed/target_info/AH_S1_target_info.txt"
df2 = read.csv(file, sep = "\t", )

df2 <- separate(df2, seqName.start.end, into = c("v1", "v2", "v3"))
df2Pos <- unite(df2, v1 , v2, col = "Pos", sep = ".") %>% select(Pos)
df2PosEnd <- unite(df2, v1 , v3, col = "Pos", sep = ".") %>% select(Pos)
df2 <- rbind(df2Pos, df2PosEnd)

df2$Pos <- as.numeric(df2$Pos)
df2 <- df2 %>% arrange(Pos) 	

df2$value <- rownames(df2)
df2$value <- as.numeric(df2$value)
df2 <- df2 %>% select(Pos, value)
df2$Gene.name <- ""
df2$Region <- "target"

df3 <- rbind(df, df2)
df3 <- df3 %>% arrange(df3$Pos) 
df3$Pos <- as.character(df3$Pos)

df3$number <- rownames(df3)

library(ggrepel)
df3 %>% 
	ggplot(aes(x=reorder(Pos, number), y=value)) +
	geom_bar(stat="identity", aes(color=Region)) +
	xlab("GRCh37 position") +
	ylab("") +
	geom_text(aes(label=Gene.name, x = Pos, y=value), hjust=-0, angle=90) +
	theme_bw()  +
	facet_grid(Region~., scales="free" )

geom_label_repel(mapping = aes(label = Gene.name, angle=45))
theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

ggsave("../output/annotation_datasets_category.pdf", width = 16, height = 8, units = "cm")
ggsave("../output/annotation_datasets_category.png", width = 16, height = 8, units = "cm")

