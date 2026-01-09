#### Jinxin Meng, jinxmeng@zju.edu.cn, 20230610, 20250309, v.2 ####
# 20231101: update function: get_freq
# 20250223: undate functions with new grammar.
# 20250223: add function: get_text_color to decide the text 'black' or 'white'.
# 20250309: add function: write_xlsx_with_comment()

library(tidyverse)

#### get_freq 筛选一个向量中出现元素次数最多的N个元素 get_element ####
# vec: 包含所有要统计元素的向量
# 过滤标准至多只能是一个参数，不设置参数输出全部结果不进行任何过滤，过滤参数需要输入数值型的值
# 默认删除向量的NA值
# 如果想要输出一个计数表格，指定out_df为T
# 如果输出一个列表，默认输出是字符串。可选择fct和num
get_freq <- function(vec, top_n = NULL, tail_n = NULL, gt = NULL, 
                        gtn = NULL, eq = NULL, ne = NULL, lt = NULL, 
                        ltn = NULL, na.rm = T, out_df = F, out_fmt = "chr") {
  args <- c(top_n, tail_n, gt, gtn, eq, ne, lt, ltn)
  if (length(args) > 1) stop("parameter error.")
  if (isTRUE(na.rm)) vec <- vec[!is.na(vec)]
  dat <- data.frame(table(vec)) %>% arrange(desc(Freq))
  if (!is.null(top_n) & is.numeric(top_n)) { out <- head(dat, n = top_n) } 
  else if (!is.null(tail_n) & is.numeric(tail_n)) { out <- tail(dat, n = top_n) }
  else if (!is.null(gt) & is.numeric(gt)) { out <- filter(dat, Freq > gt) }
  else if (!is.null(gtn) & is.numeric(gtn)) { out <- filter(dat, Freq >= gtn) }
  else if (!is.null(eq) & is.numeric(eq)) { out <- filter(dat, Freq == eq) }
  else if (!is.null(ne) & is.numeric(ne)) { out <- filter(dat, Freq != ne) }
  else if (!is.null(lt) & is.numeric(lt)) { out <- filter(dat, Freq < lt) }
  else if (!is.null(ltn) & is.numeric(ltn)) { out <- filter(dat, Freq <= ltn) }
  else { out <- dat }
  out <- dplyr::rename(out, name = vec, freq = Freq) %>% mutate(name = as.character(name))
  if (isFALSE(out_df)) {
    out <- unlist(dplyr::select(out, name), use.names = F)
    if (out_fmt == "chr") { out <- as.character(out) }
    else if (out_fmt == "num") {out <- as.numeric(out) }
  } 
  return(out)
}

#### 多位数向上或者向下取整 ####
floor_n <- function(x, n = 2) { 
  y <- x - as.numeric(str_sub(as.character(x), start = -n))
  return(y)
}

ceiling_n <- function(x, n = 2) {
  y <- x + 10^n - as.numeric((str_sub(as.character(x), start = -n)))
  return(y)
}

#### intersect_multiple 多个数据集的交集 ####
# 输入一个list，list下每个元素是一个向量
# 可使用purrr::reduce(list, \(x, y) intersect(x, y))
intersect_multiple <- function(x) {
  return(purrr::reduce(x, \(x, y) intersect(x, y)))
}

#### stat_vec ####
tools_stat_vec <- function(x, top_n = NULL, other_name = "Other", out_vec = F) {
  data <- as.character(x) %>% 
    table() %>% 
    as.data.frame(stringsAsFactors = F) %>% 
    dplyr::rename(name = 1, n = 2)
  
  if (!is.null(top_n) & is.numeric(top_n)) {
    vec <- data %>% arrange(desc(n)) %>% 
      head(n = top_n - 1) %>% 
      select(name) %>% 
      unlist(use.names = F)
    data <- data %>% 
      mutate(name = ifelse(name %in% vec, name, other_name)) %>% 
      group_by(name) %>% summarise(n = sum(n)) %>% 
      arrange(desc(n))
  }
  if (isTRUE(out_vec)) data <- data$name
  return(data)
}

#### get_text_color ####
get_text_color <- function(color) {
  rgb <- col2rgb(color) # 将颜色转换为 RGB 值
  # 计算亮度（Luminance）公式：0.299*R + 0.587*G + 0.114*B
  luminance <- apply(rgb, 2, \(x) (0.299 * x[1] + 0.587 * x[2] + 0.114 * x[3]) / 255)
  out <- purrr::map_vec(luminance, \(x) ifelse(x >= .5, 'black', 'white'))
}

#### calcu_adjusted_r2 ####
calcu_adjusted_r2 <- function(adonis_object) {
  n_observations <- adonis_object$Df[3]+1
  d_freedom <- adonis_object$Df[1]
  r2 <- adonis_object$R2[1]
  adjusted_r2 <- vegan::RsquareAdj(r2, n_observations, d_freedom)
  return(adjusted_r2)
}

#### write_xlsx_with_comment ####
write_xlsx_with_comment <- function(data, filename, comment = '###') {
  library(openxlsx)
  redStyle <- createStyle(fontColour = "red")
  
  if (is.character(comment)) {
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet1")
    writeData(wb, sheet = "Sheet1", x = comment, startRow = 1, startCol = 1)
    addStyle(wb, sheet = "Sheet1", style = redStyle, rows = 1, cols = 1)
    writeData(wb, sheet = 'Sheet1', x = data, startRow = 3)
    saveWorkbook(wb, file = filename, overwrite = T)
    
  } else if (is.list(comment)) {
    len <- length(comment)
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet1")
    
    for (i in 1:len) {
      writeData(wb, sheet = "Sheet1", x = comment[[i]], startRow = i, startCol = 1)
      addStyle(wb, sheet = "Sheet1", style = redStyle, rows = i, cols = 1)
    }
    
    writeData(wb, sheet = 'Sheet1', x = data, startRow = len + 2)
    saveWorkbook(wb, file = filename, overwrite = T)
  }
}

#### read_xlsx_multiple ####
read_xlsx_multiple <- function(file) {
  .sheet_name <- openxlsx::getSheetNames(file)
  data <- purrr::map(.sheet_name, ~ openxlsx::read.xlsx(file, sheet = .x)) %>% 
    set_names(.sheet_name)
  return(data)
}

#### ggplot_theme ####
.theme_bw_clean <- theme(axis.ticks = element_line(linewidth = .5),
                         axis.ticks.length = unit(2, 'mm'),
                         axis.text = element_text(size = 12, color = 'black'),
                         axis.title = element_text(size = 12, color = 'black'),
                         strip.background = element_blank(), 
                         plot.title = element_text(size = 13, color = 'black', face = 'bold', hjust = .5), 
                         plot.background = element_blank(),
                         panel.grid.major = element_line(linewidth = .5),
                         panel.grid.minor = element_blank(),
                         panel.border = element_rect(fill = 'NA', linewidth = .5, color = 'black'),
                         panel.background = element_blank(), 
                         panel.spacing = unit(0, 'mm'))
