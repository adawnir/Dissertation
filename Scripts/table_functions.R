ContinuousTest3=function(x, y){
  # By status
  out=NULL
  for (k in c("LUX","FRA","GS")){
    mymean=formatC(10^mean(x[y==k], na.rm=TRUE), format="f", digits=2)
    mysd=formatC(10^sd(x[y==k], na.rm=TRUE), format="f", digits=2)
    out=c(out, paste0(mymean, " (", mysd, ")"))
  }
  
  # one-way ANOVA
  mytest=aov(x ~ y, test = 'Chisq')
  mypval=formatC(mytest$`Pr(>Chi)`[2], format="e", digits=2)
  out=c(out, mypval)
  
  return(out)
}

ContinuousTest2=function(x, y){
  # By status
  out=NULL
  for (k in c("LUX","GS")){
    mymean=formatC(10^mean(x[y==k], na.rm=TRUE), format="f", digits=2)
    mysd=formatC(10^sd(x[y==k], na.rm=TRUE), format="f", digits=2)
    out=c(out, paste0(mymean, " (", mysd, ")"))
  }
  
  # T-test
  mytest=t.test(x=x[which(y=="LUX")], y=x[which(y=="GS")], alternative="two.sided")
  mypval=formatC(mytest$p.value, format="e", digits=2)
  out=c(out, mypval)
  
  return(out)
}

ReformatScientificNotation=function(x){
  xsub=x[-which(is.na(x)|!grepl("e",x))]
  mysplit=strsplit(xsub, split="e")
  xsub=sapply(mysplit, FUN=function(tmp){paste0(tmp[1],"x10_[",as.numeric(tmp[2]),"]")})
  x[-which(is.na(x)|!grepl("e",x))]=xsub
  return(x)
}

SaveExcelWithSuperscripts=function(dt, filename){
  wb <- openxlsx::createWorkbook() # create workbook
  openxlsx::addWorksheet(wb, sheetName = "data") # add sheet
  openxlsx::writeData(wb, sheet=1, x=dt, xy=c(1, 1)) # write data on workbook
  
  for(i in grep("\\_\\[([A-z0-9\\s\\-]*)\\]", wb$sharedStrings)){
    # if empty string in superscript notation, then just remove the superscript notation
    if(grepl("\\_\\[\\]", wb$sharedStrings[[i]])){
      wb$sharedStrings[[i]] <- gsub("\\_\\[\\]", "", wb$sharedStrings[[i]])
      next # skip to next iteration
    }
    
    # insert additioanl formating in shared string
    wb$sharedStrings[[i]] <- gsub("<si>", "<si><r>", gsub("</si>", "</r></si>", wb$sharedStrings[[i]]))
    
    # find the "_[...]" pattern, remove brackets and udnerline and enclose the text with superscript format
    wb$sharedStrings[[i]] <- gsub("\\_\\[([A-z0-9\\s\\-]*)\\]", "</t></r><r><rPr><vertAlign val=\"superscript\"/></rPr><t xml:space=\"preserve\">\\1</t></r><r><t xml:space=\"preserve\">", wb$sharedStrings[[i]])
  }
  openxlsx::saveWorkbook(wb, file=filename, overwrite = TRUE)
}