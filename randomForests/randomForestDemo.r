# Nicholas Ducharme-Barth
# Zach Siders
# April 3, 2018
# Brief machine learning tutorial: random forests!

# define helper function that takes true and predicted values and calculates performance metrics
	rf.performance.metrics = function(true,pred)
	{
		val = data.frame(true=true, pred=pred)
		# resulting confusion matrix
			cmx = matrix(c(sum(true==0 & pred==0),
			               sum(true==0 & pred==1),
			               sum(true==1 & pred==0),
			               sum(true==1 & pred==1)), 2, 2)
			colnames(cmx) = c("T.0","T.1")
			rownames(cmx) = c("P.0","P.1")
		# performance metrics
			PPV = cmx[2,2]/sum(cmx[2,]) # Positive Predictive Value
			NPV = cmx[1,1]/sum(cmx[1,]) # Negative Predictive Value
			sensitivity = cmx[2,2]/sum(cmx[,2])
			specificity = cmx[1,1]/sum(cmx[,1])
		return(list(cmx=cmx, PPV=PPV,NPV=NPV,sensitivity=sensitivity,specificity=specificity))	
	}

# load/install.packages
	library(randomForest)
	# library(PresenceAbsence)

# bring in the data
# we will use the train and test datasets from the Kaggle machine learning competition
# these can be found at
# https://www.kaggle.com/c/titanic/data#

	# path.read = "/your_file_path/"
	# path.read = "/Users/nicholasducharme-barth/Desktop/Course Work/Spring18/RandomForestDemo/"
	path.read = "/Volumes/HDD/Users/Zach/Documents/AAA_Git_projects/uf_lab/randomForests/"
	train = read.csv(paste0(path.read,"train.csv"))
	test = read.csv(paste0(path.read,"test.csv"))

# look at our data
	str(train)
	str(test)

# subset our data to desired columns and remove NAs
	train.df = na.omit(train[,c(2,3,5,6,7,8,10)])
	test.df = na.omit(test[,c(2,4,5,6,7,9)])

# run a random forest model
# convert to factor for classification
	set.seed(7890)
	simple.rf = randomForest(as.factor(Survived) ~ Sex, data=train.df)
	complicated.rf = randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df)

# get predictions
	simple.pred = predict(simple.rf, test.df)
	complicated.pred = predict(complicated.rf, test.df)

# look at and evaluate output...
	partialPlot(simple.rf, train.df, "Sex",0)

	# png(file=paste0(path.read,"Response.curves.png"), width=600, height=400)
		par(mfrow=c(2,3), oma=c(0.4,2,0.4,0.4))
		partialPlot(complicated.rf, train.df, "Sex",0)
		mtext("Perished", side=2, line=3)
		partialPlot(complicated.rf, train.df, "Age",0)
		partialPlot(complicated.rf, train.df, "Fare",0)
		partialPlot(complicated.rf, train.df, "Sex",1)
		mtext("Survived", side=2, line=3)
		partialPlot(complicated.rf, train.df, "Age",1)
		partialPlot(complicated.rf, train.df, "Fare",1)
	# dev.off()
# how does increasing trees impact performance
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,ntree=10)
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,ntree=100)
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,ntree=500)
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,ntree=1000)

# variables at each split?
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,mtry=3)
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,mtry=4)
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,mtry=5)
	randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train.df,mtry=6)

# variable importance
	varImpPlot(simple.rf, main="Variable Importance for simple")

	varImpPlot(complicated.rf, main="Variable Importance for complicated")

# need to make a new train/test dataset
	sample = sample(1:nrow(train.df),floor(nrow(train.df)*0.8), replace=FALSE)
	train2 = train.df[sample,]
	test2 = train.df[-sample,]

#rerun RFs with subsetted datasets
	simple.rf2 = randomForest(as.factor(Survived) ~ Sex, data=train2)
	complicated.rf2 = randomForest(as.factor(Survived) ~ Pclass + Sex + Age + SibSp + Parch + Fare, data=train2)

# generate predictions
	simple.pred2 = predict(simple.rf2, test2)
	complicated.pred2 = predict(complicated.rf2, test2)

#looking at the votes
	(complicated.pred.votes = predict(complicated.rf2, test2, type='prob'))

# Performance metrics
	rf.performance.metrics(test2$Survived, simple.pred2)
	rf.performance.metrics(test2$Survived, complicated.pred2)

# Receiver Operator Characteristic Curve
	rocr <- function(forest,samples,name){
		library(ROCR)
		cls_pred <- prediction(predict(forest, type='prob',samples)[,2], samples$Survived)
		cls_perf <- performance(cls_pred, "tpr", "fpr")
		cls_nperf <- performance(cls_pred, "tnr", "fnr")
		cls_ss <- performance(cls_pred, "sens", "spec")
		cls_ss@x.values[[1]] <- 1-cls_ss@x.values[[1]]
		cls_auc <- performance(cls_pred, measure='auc')@y.values[[1]]
		# plot(cls_perf, main="ROC Curve for Random Forest", col=2, lwd=2)
		# abline(a=0, b=1, lwd=2, lty=2, col="gray")
			plot(cls_ss, main=paste0("ROC Curve for ",name), col=2, lwd=2, xlab="1-Specificity")
			abline(a=0, b=1, lwd=2, lty=2, col="gray")
			text(0.8,0.2,paste0("AUC = ", round(cls_auc,3)), cex=1.2)

		return(list(tpr=cls_perf, ss=cls_ss, tnr=cls_nperf, auc=cls_auc))
	}

# ROC for test data
	roc.simp.test <- rocr(simple.rf2, test2, "simple")
	roc.comp.test <- rocr(complicated.rf2, test2, "complex")

# ROC for train data
	roc.simp.train <- rocr(simple.rf2, train2, "simple")
	roc.comp.train <- rocr(complicated.rf2, train2, "complex")

#can just use the original train DF since we subsetted
	roc.simp.both <- rocr(simple.rf2, train, "simple")
	roc.comp.both <- rocr(complicated.rf2, train, "complex")


# Threshold of Maximum Sensitivity - Specificity
	thres <- roc.comp.test$ss@alpha.values[[1]][which.max(((1-roc.comp.test$ss@x.values[[1]]) + roc.comp.test$ss@y.values[[1]]) - 1)]
	thres_mat <- which.max(((1-roc.comp.test$ss@x.values[[1]]) + roc.comp.test$ss@y.values[[1]]) - 1)

# plot up the complicated train/test ROC
	# png(file=paste0(path.read,"ROC.curve.comp.png"), width=400, height=400)
		plot(roc.comp.test$ss, main="ROC Curve for complicated", col="dodgerblue", lwd=2, xlab="1-Specificity")
		lines(roc.comp.train$ss@x.values[[1]],roc.comp.train$ss@y.values[[1]], col='seagreen4', lwd=2)
		points(roc.comp.test$ss@x.values[[1]][thres_mat], roc.comp.test$ss@y.values[[1]][thres_mat], col='black', pch=4, cex=1.3)
		abline(a=0, b=1, lwd=2, lty=2, col="gray")
		text(0.8,0.2,paste0("Test AUC = ", round(roc.comp.test$auc,3)), cex=1.2, col="dodgerblue")
		text(0.8,0.3,paste0("Train AUC = ", round(roc.comp.train$auc,3)), cex=1.2, col="seagreen4")
	# dev.off()




