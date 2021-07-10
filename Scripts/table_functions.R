ContinuousTest3=function(x, y){
  # By status
  out=NULL
  for (k in c("LUX","FRA","GS")){
    mymean=formatC(mean(x[y==k], na.rm=TRUE), format="f", digits=2)
    mysd=formatC(sd(x[y==k], na.rm=TRUE), format="f", digits=2)
    out=c(out, paste0(mymean, " (", mysd, ")"))
  }
  
  # T-test
  mytest1=t.test(x=x[which(y=="LUX")], y=x[which(y=="FRA")], alternative="two.sided")
  mytest2=t.test(x=x[which(y=="LUX")], y=x[which(y=="GS")], alternative="two.sided")
  mytest3=t.test(x=x[which(y=="FRA")], y=x[which(y=="GS")], alternative="two.sided")
  mypval1=formatC(mytest1$p.value, format="e", digits=2)
  mypval2=formatC(mytest2$p.value, format="e", digits=2)
  mypval3=formatC(mytest3$p.value, format="e", digits=2)
  out=c(out, mypval1, mypval2, mypval3)
  
  return(out)
}

ContinuousTest2=function(x, y){
  # By status
  out=NULL
  for (k in c("LUX","GS")){
    mymean=formatC(mean(x[y==k], na.rm=TRUE), format="f", digits=2)
    mysd=formatC(sd(x[y==k], na.rm=TRUE), format="f", digits=2)
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
  wb <- createWorkbook() # create workbook
  addWorksheet(wb, sheetName = "data") # add sheet
  writeData(wb, sheet=1, x=dt, xy=c(1, 1)) # write data on workbook
  
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
  saveWorkbook(wb, file=filename, overwrite = TRUE)
}