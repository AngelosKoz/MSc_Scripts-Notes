####### EXERCISE 2 #######
# # # # IMPORTANT NOTICE: paths need to be changed to corresponds to file locations # # # #
# If Assignment2_Kozwnakis_Angelos is in the same folder, this will work, else it has to be specified with sys.path.append which has been commented out

import pandas as pd
import numpy as np
# Load the Naive Bayes Classifier train and predict
# Since my default directory is set differently, i had to explicitly specify the path to load the function
# To avoid having an even larger scripts, i chose to import the functions instead of adding them. 
# sys append needs to be changed to appropriate path if previous .py file is not in the same folder
#import sys
#sys.path.append("/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/")
from Assignment2_Kozwnakis_Angelos import train_NBC, predict_NBC
# Import logistic regression, one hot encoder and random spliting from sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
import matplotlib.pyplot as plt

# Create the trivial train function
def trivial_train(Y):
    #from collections import Counter
    # We find the unique classes with set, and then count them, keeping the maximum counter number. (Most frequent)
    #print(Counter(Y)) # We can use this to see the occurences of each class. Since it is not implied that we can use external libraries, i just used it for reference
    return float(max(set(list(Y)), key=list(Y).count)) #Since .count does not work for numpy, we use the list

# Create the trivial predict function.
def trivial_predict(X, trivial_class):
    # All this function will do, is assign all the samples to the most common/frequent class.
    # We can do this by creating an array with length equal to the samples and fill it with the class from trivial train
    return np.full(shape = (len(X),), fill_value=trivial_class)

def run_analysis(Xcat, Ycat, Xcont, Ycont, Xmix = None, Ymix = None, L=1, num_runs=100, size=0.25):
    
    #Initialize all the keys for our different data types and classifiers
    accuracies = {}
    accuracies['continuous'] = {}
    accuracies['continuous']['triv'] = {}
    accuracies['continuous']['NBC'] = {}
    accuracies['continuous']['LR'] = {}
    accuracies['categorical'] = {}
    accuracies['categorical']['triv'] = {}
    accuracies['categorical']['triv_OHE'] = {}
    accuracies['categorical']['NBC'] = {}
    accuracies['categorical']['NBC_OHE'] = {}
    accuracies['categorical']['LR'] = {}
    accuracies['categorical']['LR_OHE'] = {}
    
    # Create dictionary values for each K for all they keys above (different classifiers and data type)
    K = list(range(50, 101, 10))
    for k in K:
            accuracies['continuous']['triv'][k] = []
            accuracies['continuous']['NBC'][k] = []
            accuracies['continuous']['LR'][k] = []
            accuracies['categorical']['triv'][k] = []
            accuracies['categorical']['triv_OHE'][k] = []
            accuracies['categorical']['NBC'][k] = []
            accuracies['categorical']['NBC_OHE'][k] = []
            accuracies['categorical']['LR'][k] = []
            accuracies['categorical']['LR_OHE'][k] = []
    
    #### PRE-PROCESSING ####
    # Create categorical data without One Hot Encoding (OHE)
    D_cat = (np.array([len(np.unique(Xcat[:,i])) for i in range(Xcat.shape[1])])).ravel() 
    OHE = OneHotEncoder()
    # Apply OHE to the categorical X data (the features of Dataset A)
    Xcat_OHE = OHE.fit_transform(Xcat).toarray()
    # Extracting the dimensions of the one-hot encoded dataset (mainly the number of distinct values after encoding)
    D_cat_OHE = (np.array([int(np.max(Xcat_OHE[:,i])+1) for i in range(Xcat_OHE.shape[1])])).ravel()
    
    # Start iterations
    for _ in range(num_runs):
        # Split the data into (75% training set and 25% test set)
        # We do it for categorical and continuous seperately:
        
        # 1) For Categorical:
            # a) Without OHE
        X_train_cat, X_test_cat, Y_train_cat, Y_test_cat = train_test_split(Xcat, Ycat, test_size=size)
            # b) With OHE
        X_train_cat_OHE, X_test_cat_OHE, Y_train_cat_OHE, Y_test_cat_OHE = train_test_split(Xcat_OHE, Ycat, test_size=size)
        # 2) For continuous
        X_train_cont, X_test_cont, Y_train_cont, Y_test_cont = train_test_split(Xcont, Ycont, test_size=size)
    
        # Now we start iterations for different % of starting training set kept.
        for k in K:
            #### 1) CATEGORICAL:
                        
            ##  A) Without OHE
            X_train_cat_K = X_train_cat[: int( (k/100) * len(X_train_cat))]
            Y_train_cat_K = Y_train_cat[: int( (k/100) * len(Y_train_cat))]
            
                # Trivial classifier
            predict_cat_triv = trivial_predict(X_test_cat, trivial_train(Y_train_cat_K))
            correct_predictions_cat_triv = np.sum(predict_cat_triv == Y_test_cat)
            total_predictions_cat_triv = len(Y_test_cat)            
            accuracy_cat_triv = correct_predictions_cat_triv/total_predictions_cat_triv
            accuracies['categorical']['triv'][k].append(accuracy_cat_triv)
            #print(f'Samples: {len(X_train_cat_K)}')
            #print(f'Correct pred: {correct_predictions_cat_triv}')
            #print(f'Total pred: {total_predictions_cat_triv}')
            #print(f'Accuracy: {accuracy_cat_triv}\n')
                # Naive Bayes Classifier
            model_cat_NBC = train_NBC(X = X_train_cat_K, X_dtype="categorical", Y = Y_train_cat_K, L = L, D_categorical=D_cat)
            predict_cat_NBC, classprobs_cat = predict_NBC(model = model_cat_NBC, X = X_test_cat, X_dtype="categorical")
            correct_predictions_cat_NBC = np.sum(predict_cat_NBC == Y_test_cat)
            total_predictions_cat_NBC = len(Y_test_cat)
            accuracy_cat_NBC = correct_predictions_cat_NBC/total_predictions_cat_NBC
            accuracies['categorical']['NBC'][k].append(accuracy_cat_NBC)
                # Logistic Regression
            model_cat_LR = LogisticRegression(C = 1/L, max_iter = 10000)
            model_cat_LR.fit(X_train_cat_K, Y_train_cat_K)
            predict_cat_LR = model_cat_LR.predict(X_test_cat)
            correct_predictions_cat_LR = np.sum(predict_cat_LR == Y_test_cat)
            total_predictions_cat_LR = len(Y_test_cat)
            accuracy_cat_LR = correct_predictions_cat_LR/total_predictions_cat_LR
            accuracies['categorical']['LR'][k].append(accuracy_cat_LR)
            
            ##  B) With OHE
            X_train_cat_K_OHE = X_train_cat_OHE[: int( (k/100) * len(X_train_cat_OHE))]
            Y_train_cat_K_OHE = Y_train_cat_OHE[: int( (k/100) * len(Y_train_cat_OHE))]
            
                # Trivial classifier
            predict_cat_triv_OHE = trivial_predict(X_test_cat_OHE, trivial_train(Y_train_cat_OHE))
            correct_predictions_cat_triv_OHE = np.sum(predict_cat_triv_OHE == Y_test_cat_OHE)
            total_predictions_cat_triv_OHE = len(Y_test_cat_OHE)
            
            accuracy_cat_triv_OHE = correct_predictions_cat_triv_OHE/total_predictions_cat_triv_OHE
            accuracies['categorical']['triv_OHE'][k].append(accuracy_cat_triv_OHE)
                # Naive Bayes Classifier
            model_cat_NBC_OHE = train_NBC(X = X_train_cat_K_OHE, X_dtype="categorical", Y = Y_train_cat_K_OHE, L = L, D_categorical=D_cat_OHE)
            predict_cat_NBC_OHE, classprobs_cat_OHE = predict_NBC(model = model_cat_NBC_OHE, X = X_test_cat_OHE, X_dtype="categorical")
            correct_predictions_cat_NBC_OHE = np.sum(predict_cat_NBC_OHE == Y_test_cat_OHE)
            total_predictions_cat_NBC_OHE = len(Y_test_cat_OHE)
            accuracy_cat_NBC_OHE = correct_predictions_cat_NBC_OHE/total_predictions_cat_NBC_OHE
            accuracies['categorical']['NBC_OHE'][k].append(accuracy_cat_NBC_OHE)
                # Logistic Regression
            model_cat_LR_OHE = LogisticRegression(C = 1/L, max_iter = 10000)    
            model_cat_LR_OHE.fit(X_train_cat_K_OHE, Y_train_cat_K_OHE)
            predict_cat_LR_OHE = model_cat_LR_OHE.predict(X_test_cat_OHE)
            correct_predictions_cat_LR_OHE = np.sum(predict_cat_LR_OHE == Y_test_cat_OHE)
            total_predictions_cat_LR_OHE = len(Y_test_cat_OHE)
            accuracy_cat_LR_OHE = correct_predictions_cat_LR_OHE/total_predictions_cat_LR_OHE
            accuracies['categorical']['LR_OHE'][k].append(accuracy_cat_LR_OHE)

            ### 2) CONTINUOUS
            X_train_cont_K = X_train_cont[: int( (k/100) * len(X_train_cont))]
            Y_train_cont_K = Y_train_cont[: int( (k/100) * len(Y_train_cont))]
            
                # Trivial classifier
            predict_cont_triv = trivial_predict(X_test_cont, trivial_train(Y_train_cont_K))
            correct_predictions_cont_triv = np.sum(predict_cont_triv == Y_test_cont)
            total_predictions_cont_triv = len(Y_test_cont)
            accuracy_cont_triv = correct_predictions_cont_triv/total_predictions_cont_triv
            accuracies['continuous']['triv'][k].append(accuracy_cont_triv)
                # Naive Bayes Classifier
            model_cont_NBC = train_NBC(X = X_train_cont_K, X_dtype="continuous", Y = Y_train_cont_K)
            predict_cont_NBC, classprobs_cont = predict_NBC(model = model_cont_NBC, X = X_test_cont, X_dtype="continuous")
            correct_predictions_cont_NBC = np.sum(predict_cont_NBC == Y_test_cont)
            total_predictions_cont_NBC = len(Y_test_cont)
            accuracy_cont_NBC = correct_predictions_cont_NBC/total_predictions_cont_NBC
            accuracies['continuous']['NBC'][k].append(accuracy_cont_NBC)
                # Logistic Regression
            model_cont_LR = LogisticRegression(C = 1/L, max_iter = 10000)
            model_cont_LR.fit(X_train_cont_K, Y_train_cont_K)
            predict_cont_LR = model_cont_LR.predict(X_test_cont)
            correct_predictions_cont_LR = np.sum(predict_cont_LR == Y_test_cont)
            total_predictions_cont_LR = len(Y_test_cont)
            accuracy_cont_LR = correct_predictions_cont_LR/total_predictions_cont_LR
            accuracies['continuous']['LR'][k].append(accuracy_cont_LR)
          
    #mean_accuracies = {k1: {k2: {k3: np.mean(v3) for k3, v3 in v2.items()} for k2, v2 in v1.items()} for k1, v1 in accuracies.items()}
    # Calculate mean accuracy for all data types and classifiers
    mean_accuracies = {
        k1: {
            k2: {
                k3: np.mean(v3) for k3, v3 in v2.items()
            } for k2, v2 in v1.items()
        } for k1, v1 in accuracies.items()
    }

    #Create a combined plot showing the relationship between accuracy and sample size for each classifier
    # First we access the different data types and their classifiers
    for dtype, classifier in mean_acc.items():
        # Create a new figure for each dataset
        plt.figure(figsize=(10, 6))  
        # Iterate through the classifiers and the different training set size
        for classifier, sample_sizes in classifier.items():
            # Extract sample sizes and their corresponding accuracies
            sample_sort = sorted(sample_sizes.items())  # Sort by sample size
            train_percent = [k for k, v in sample_sort] #Keep the sample sizes
            accuracies = [v for k, v in sample_sort] #Keep mean accuracy values

            plt.plot(train_percent, accuracies, label=classifier)
        
        plt.title(f'Accuracy vs Training Sample Size  ({dtype})')
        plt.xlabel('Sample Size of Training Set (%)')
        plt.ylabel(f'Average Accuracy (Reps: {num_runs})')
        plt.legend()
        plt.grid(True)
        plt.xticks(train_percent)
        plt.show()


    return mean_accuracies, accuracies


X_cat = (pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/Assignment3_DataCode/Dataset3.2_A_X.csv', header = None)).values
Y_cat = ((pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/Assignment3_DataCode/Dataset3.2_A_Y.csv', header = None)).values).ravel()
X_cont = (pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/Assignment3_DataCode/Dataset3.2_B_X.csv', header = None, delimiter=";")).values
Y_cont = ((pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/Assignment3_DataCode/Dataset3.2_B_Y.csv', header = None)).values).ravel()

mean_acc, acc = run_analysis(Xcat = X_cat, Ycat = Y_cat, Xcont = X_cont, Ycont = Y_cont, L=1, num_runs=100)

# We can load the dictionary like so:
#import ast
#with open('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/dict_mean_1k', 'r') as f:
#    acc_dict = f.read()
#    mean_accuracies = ast.literal_eval(acc_dict)



####### EXERCISE 3 #######
X_3_3 = (pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/Assignment3_DataCode/Dataset3.3_X.csv', header = None)).values
Y_3_3 = ((pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_3/Assignment3_DataCode/Dataset3.3_Y.csv', header = None)).values).ravel()

# L1 regularization cannot use 'lbfgs' (default for LR) because it does not support L1 penalty. 
# So we use 'saga' and keep it for the no penalty as well for consistency

# No penalization
LR_no_pen = LogisticRegression(penalty =None, solver = 'saga', max_iter = 10000)
LR_no_pen.fit(X_3_3, Y_3_3)

# C = 1 / λ
# Lasso regularization: λ = 0.5 (C = 2)
LR_lasso_0_5 = LogisticRegression(penalty='l1', C=2, solver='saga', max_iter=10000)
LR_lasso_0_5.fit(X_3_3, Y_3_3)
# Lasso regularization: λ = 10 (C = 0.1)
LR_lasso_10 = LogisticRegression(penalty='l1', C = 1/10, solver='saga', max_iter=10000)
LR_lasso_10.fit(X_3_3, Y_3_3)
# Lasso regularization: λ = 100 (C = 0.01)
LR_lasso_100 = LogisticRegression(penalty='l1', C = 1/100, solver='saga', max_iter=10000)
LR_lasso_100.fit(X_3_3, Y_3_3)

# Weights of each model
# Each weight array here has a length equal to the features (in our case we have 1000x30 initial matrix, so 30 features)
# The absolute value of the weight shows how much the feature contributes to the model's predictions (higher absolute values, higher impact)
weights_no_penalty = LR_no_pen.coef_[0]
weights_lasso_0_5 = LR_lasso_0_5.coef_[0]
weights_lasso_10 = LR_lasso_10.coef_[0]
weights_lasso_100 = LR_lasso_100.coef_[0]

# Combine the weights
all_weights = np.hstack([weights_no_penalty, weights_lasso_0_5, weights_lasso_10, weights_lasso_100])

# Later on the plots did not have negative values, so this had to be implemented to have the same plots. Better for comparison
# Basicly we manually apply limits for y_axis
ymin, ymax = all_weights.min(), all_weights.max()
yrange = ymax - ymin
ymin -= 0.1 * yrange  # Add some space at the bottom
ymax += 0.1 * yrange  # Add some space at the top

# We can combine the models and use a loop for the plots
models = [
    ('No Penalty', weights_no_penalty),
    ('Lasso λ=0.5', weights_lasso_0_5),
    ('Lasso λ=10', weights_lasso_10),
    ('Lasso λ=100', weights_lasso_100)
]

for name, weights in models:
    print(weights)
    plt.figure(figsize=(10, 5))
    plt.bar(range(len(weights)), weights)
    plt.ylim(ymin, ymax)  # Same y-axis for comparison
    plt.xlabel('Weight Index')
    plt.ylabel('Weight Value')
    plt.title(f'Logistic Regression Weights ({name})') 
    plt.show()
