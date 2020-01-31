#' Visualized network
#' @import shiny
#' @importFrom data.table melt
#' @importFrom grDevices pdf
#' @importFrom stats na.omit
#' @param Auto_WGCNA_Output R object with GMIC bayesian network
#' @param Filter_unconnected_ME a logical value. If TRUE, the default,
#' unconnected modules will be removed from the final network. If FALSE,
#' all modules will be shown.
#' @return a shiny object for network visualization. 
#' @examples GMIC_Final_dir<-system.file("extdata", "GMIC_Final.Rdata", 
#' package = "GmicR", mustWork = TRUE)
#' load(GMIC_Final_dir)
#' if(interactive()){
#' Gmic_viz(GMIC_Final)}
#' @export

Gmic_viz<-function(Auto_WGCNA_Output,
Filter_unconnected_ME = TRUE){

# setting up tables    
Gmic_Output<-Auto_WGCNA_Output$Output_tabu_net
Module_Annotations<-Auto_WGCNA_Output$GO_Query
Summary_Table<-Auto_WGCNA_Output$GO_table
GSEAGO_query<-Auto_WGCNA_Output$GSEAGO_Builder_Output$GSEAGO

# Module_Annotations subsetting frame
Module_selection_frame<-Summary_Table
Module_selection_frame$Name<-paste(rownames(Module_selection_frame), 
Module_selection_frame$GO_module_name, sep = ": ")
Module_selection_frame$GSEAGO_ID<-c(Module_selection_frame$modules+1)

# annotating GSEAGO_query table
table_GSEAGO_query<-data.table::melt(GSEAGO_query,
id.vars=c("GOBPID", "Pvalue", "OddsRatio", 
"ExpCount", "Count", "Size", 
"Term","Genes"))
table_GSEAGO_query$Module<-c(table_GSEAGO_query$L1-1)


# preparing for highlighting module name
table_GSEAGO_query$Module_Name<-c("")
MSF_ids<-rownames(Module_selection_frame)
for(i in MSF_ids){
table_GSEAGO_query[table_GSEAGO_query$Term == Module_selection_frame[i,
]$GO_module_name & 
table_GSEAGO_query$Module == Module_selection_frame[i,
]$modules,c("Module_Name")]<-c("Yes")
}

table_GSEAGO_query<-table_GSEAGO_query[c("Module_Name",
"GOBPID", "Pvalue",
"OddsRatio", "ExpCount", "Count", "Size", 
"Term","Genes", "Module")]



# formatting Module_Annotations
if(!is.null(Module_Annotations$kTotal)){
Module_Annotations$kTotal<-round(Module_Annotations$kTotal,3)
Module_Annotations$kWithin<-round(Module_Annotations$kWithin,3)
}

# connected nodes
if(Filter_unconnected_ME==TRUE){
arc_nodes<-arcs(Gmic_Output$ave.boot)
q_nodes<-unique(c(arc_nodes[,1],arc_nodes[,2]))
}else if(Filter_unconnected_ME==FALSE){
q_nodes<-nodes(Gmic_Output$ave.boot)
}

# mapping to module table
sel_df<-Summary_Table[c(1)]
sel_df$Name<-paste(rownames(sel_df), 
sel_df$GO_module_name, sep = ": ")
sel_df$Value<-rownames(sel_df)

# Connected MEs
module_map<-sel_df[rownames(sel_df) %in% q_nodes,]
# only cell sigs influencing MEs
True_GMIC<-Batch_Net(Gmic_Output, 
Node_ids = sel_df$Value, 
relationship_type = "mb")

# getting nodes
True_GMIC_n<-nodes(True_GMIC$ave.boot)

# GMIC data
GMIC<-data.frame(Name=True_GMIC_n,
Value=True_GMIC_n, 
stringsAsFactors = FALSE)

# cell sig selection
cell_map<-data.frame(Name=True_GMIC_n[!True_GMIC_n %in% sel_df$Value],
Value=True_GMIC_n[!True_GMIC_n %in% sel_df$Value], 
stringsAsFactors = FALSE)

cell_map<-cell_map[order(cell_map$Name),]
rownames(cell_map)<-NULL

# args for radio label
if(length(cell_map$Value)!=0){
cell_option_label<-"Non-MEs variable selection:"
} else if(length(cell_map$Value)==0){
cell_option_label<-" " }

GMIC_net_query<-data.frame(Name=True_GMIC_n,
Value=True_GMIC_n, 
stringsAsFactors = FALSE)


# UI ------------------------------------------------------------------

ui <- fluidPage(
sidebarLayout(sidebarPanel(width = 3,


# tab1 GMIC -----------------------------------------------------------
conditionalPanel('input.dataset === "GMIC_net_query"',
downloadButton(outputId = "down1", label = "Download the plot"),


# Input: Checkbox for query node color  ----
radioButtons(inputId = "color",
label = "Query node color:",
inline  = TRUE,
choices =  c("lightblue", "pink", "white"),
selected =  c("lightblue","pink", "white")[1]),

# Input: Checkbox for inverse arc colors  ----
radioButtons(inputId = "edgecolor",
label = "inverse edge color selection:",
inline  = TRUE,
choices =  c("black","grey", "red"),
selected =  c("black","grey", "red")[1]), 

# Input: Checkbox for inverse arc line pattern  ----
radioButtons(inputId = "line",
label = "inverse edge pattern selection:",
inline  = TRUE,
choices =    c("longdash","blank", "solid", "dotted" ),
selected =    c("longdash","blank", "solid", "dotted" )[1]),  

# Input: Checkbox for query node relatioships ----
radioButtons(inputId = "relationship",
label = "selected node relationships:",
inline  = TRUE,
choices =  c("mb", "nbr"),
selected =  c("mb", "nbr")[1]),

# Input: Query node selection ----
actionLink("selectall","Select All"),

checkboxGroupInput(inputId = "module",
label = "GMIC module query selection:",
inline  = TRUE,
choiceValues = module_map$Value,
choiceNames = module_map$Name),

# selectallCells
actionLink("selectallCells","Select All"),

# Input: Query cel signature selection ----
checkboxGroupInput(inputId = "cell_map",
label = cell_option_label,
inline  = TRUE,
choiceValues = cell_map$Value,
choiceNames = cell_map$Name,
selected = NULL)),


# tab2 module annotations ----------------------------------------------
conditionalPanel('input.dataset === "Module_Annotations"',
# Button
checkboxGroupInput(inputId="show_vars", "Columns:",
width = 2, 
choices = names(Module_Annotations), 
selected = names(Module_Annotations)),

# Input: Checkbox for inverse arc colors  ----
# Input: Query node selection ----
actionLink("selectallModules","Select All"),
checkboxGroupInput(inputId = "module_subset",
label = "Module selection:",
inline  = TRUE,
choiceNames = Module_selection_frame$Name,
choiceValues = Module_selection_frame$modules)),


# Tab3 summary table ------------------------------------------------------
conditionalPanel('input.dataset === "Summary_Table"',
helpText(" ")),



# tab4 GSEAGO table -------------------------------------------------------
conditionalPanel('input.dataset === "table_GSEAGO_query"',
# Input: radiobuttons for GSEAGO table  ----
radioButtons(inputId = "GSEAGO_query_id",
label = "GSEAGO table for:",
inline  = TRUE,
choiceNames = Module_selection_frame$Name,
choiceValues = Module_selection_frame$modules,
selected = Module_selection_frame$modules[1]))),

# mainPanels for tabs
mainPanel(tabsetPanel(id = 'dataset',
tabPanel("GMIC_net_query", plotOutput(outputId = "GMIC_net")),
tabPanel("Module_Annotations", DT::dataTableOutput("tab2_Query")),
tabPanel("Summary_Table", DT::dataTableOutput("tab3_table")),
tabPanel("table_GSEAGO_query", DT::dataTableOutput("tab4"))
))))


# server --------------------------------------------------------------


server <- function(input, output, session) {



# updating tab1 selectall ---------------------------------------------
# GMIC module selection  
observe({if(input$selectall == 0) return(NULL) 
else if (input$selectall%%2 == 0){
updateCheckboxGroupInput(session,
"module","GMIC module query selection:",
inline  = TRUE,
choiceValues = module_map$Value,
choiceNames = module_map$Name)}
else{updateCheckboxGroupInput(session,
"module","GMIC module query selection:",inline  = TRUE,
choiceValues = module_map$Value,
choiceNames = module_map$Name,
selected = module_map$Value)}})

# cell signatures
observe({if(input$selectallCells == 0) return(NULL) 
else if (input$selectallCells%%2 == 0){
updateCheckboxGroupInput(session,
"cell_map",cell_option_label,
inline  = TRUE,
choiceValues = cell_map$Value,
choiceNames = cell_map$Name)}
else{updateCheckboxGroupInput(session,
"cell_map",cell_option_label,
inline  = TRUE,
choiceValues = cell_map$Value,
choiceNames = cell_map$Name,
selected = cell_map$Value)}}) 


# updating tab2 selectall ---------------------------------------------
observe({if(input$selectallModules == 0) return(NULL) 
else if (input$selectallModules%%2 == 0){
updateCheckboxGroupInput(session,
"module_subset","Module selection:",
inline  = TRUE,
choiceValues = Module_selection_frame$modules,
choiceNames = Module_selection_frame$Name)}
else{updateCheckboxGroupInput(session,
"module_subset","Module selection:",
inline  = TRUE,
choiceValues = Module_selection_frame$modules,
choiceNames = Module_selection_frame$Name,
selected = Module_selection_frame$modules)}}) 


# Tab2 columns to display
Module_Query2 = Module_Annotations
output$tab2_Query <- DT::renderDataTable({DT::datatable(Module_Query2[
Module_Query2$modules %in% input$module_subset,
input$show_vars, drop = FALSE], filter = 'top',
extensions = 'Buttons',
options = list(dom = 'Bfrtip', pageLength = -1,
buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))) })



# customize the length drop-down menu
output$tab3_table <- DT::renderDataTable({
DT::datatable(Summary_Table,  filter = 'top',
extensions = 'Buttons',
options = list(dom = 'Bfrtip', pageLength = -1,
buttons = c('copy', 'csv', 'excel', 'pdf', 'print')))}) 

# Tab4 Modules to display
table_GSEAGO_query2 = table_GSEAGO_query
output$tab4 <- DT::renderDataTable({DT::datatable(table_GSEAGO_query2[
table_GSEAGO_query2$Module %in% input$GSEAGO_query_id,c(seq(1,9))], filter = 'top',
extensions = 'Buttons',
options = list(dom = 'Bfrtip', pageLength = -1,
buttons = c('copy', 'csv', 'excel', 'pdf', 'print'))) })

output$GMIC_net <- renderPlot({

input$module
input$cell_map
node_id <- c(input$cell_map,input$module)

lty<-input$line #"solid"
col_line<-input$edgecolor
col_fill = input$color
lwd_set <- 1
sing_sub_net<-Batch_Net(Gmic_Output, 
Node_ids = node_id, relationship_type = input$relationship)

sing_sub_net$highlight = list(nodes=node_id,
fill = col_fill,
arcs=sing_sub_net$InverseR_Output,
lty=lty, lwd= lwd_set,
col= col_line)

graphviz.plot(sing_sub_net$ave.boot,
shape = "ellipse",
highlight= sing_sub_net$highlight,
main = "")

})

## call the plot function when downloading the image
output$down1 <- downloadHandler(
filename =  function() {
paste("net", "pdf", sep=".")
},
# content is a function with argument file. content writes the plot to the device
content = function(file) {

pdf(file) # open the pdf device
input$module
input$cell_map
node_id <- c(input$cell_map,input$module)

lty<-input$line #"solid"
col_line<-input$edgecolor
col_fill = input$color
lwd_set <- 1
sing_sub_net<-Batch_Net(Gmic_Output, 
Node_ids = node_id, relationship_type = input$relationship)

sing_sub_net$highlight = list(nodes=node_id,
fill = col_fill,
arcs=sing_sub_net$InverseR_Output,
lty=lty, lwd= lwd_set,
col= col_line)

graphviz.plot(sing_sub_net$ave.boot,
shape = "ellipse",
highlight= sing_sub_net$highlight,
main = "")
dev.off()  # turn the device off
},contentType='pdf' )



}
shinyApp(ui, server)

}
