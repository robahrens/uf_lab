# Nicholas Ducharme-Barth
# April 3, 2018
# Brief machine learning tutorial: random forests!

# define helper function that takes true and predicted values and calculates performance metrics
	rf.performance.metrics = function(true,pred)
	{
		val = data.frame(true=true, pred=pred)
		# resulting confusion matrix
			cmx = cmx(val, which.model = 1, na.rm = TRUE)
		# performance metrics
			PPV = cmx[1,1]/sum(cmx[1,]) # Positive Predictive Value
			NPV = cmx[2,2]/sum(cmx[2,]) # Negative Predictive Value
			sensitivity = cmx[1,1]/sum(cmx[,1])
			specificity = cmx[2,2]/sum(cmx[,2])
		return(list(PPV=PPV,NPV=NPV,sensitivity=specificity))	
	}

# load/install.packages
	library(randomForest)
	library(PresenceAbsence)

# bring in the data
# we will use the train and test datasets from the Kaggle machine learning competition
# these can be found at
# https://www.kaggle.com/c/titanic/data#

	path.read = "your/path/here/"
	# path.read = "/Users/nicholasducharme-barth/Desktop/Course Work/Spring18/RandomForestDemo/"
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

	par(mfrow=c(2,3))
	partialPlot(complicated.rf, train.df, "Sex",0)
	partialPlot(complicated.rf, train.df, "Age",0)
	partialPlot(complicated.rf, train.df, "Fare",0)
	partialPlot(complicated.rf, train.df, "Sex",1)
	partialPlot(complicated.rf, train.df, "Age",1)
	partialPlot(complicated.rf, train.df, "Fare",1)

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

	varImpPlot(complex.rf, main="Variable Importance for complicated")

# Receiver Operator Characteristic Curve
rocr <- function(forest,samples,name){
	library(ROCR)
	cls_pred <- prediction(predict(forest, type='prob')[,2], samples)
	cls_perf <- performance(cls_pred, "tpr", "fpr")
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
roc.simp.test <- rocr(simple.rf, test, "simple")
roc.comp.test <- rocr(complex.rf, test, "complex")

roc.simp.train <- rocr(simple.rf, train, "simple")
roc.comp.train <- rocr(complex.rf, train, "complex")

roc.simp.both <- rocr(simple.rf, test, "simple")
roc.comp.both <- rocr(complex.rf, test, "complex")






