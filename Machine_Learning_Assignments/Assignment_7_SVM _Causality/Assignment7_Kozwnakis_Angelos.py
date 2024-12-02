from sklearn.datasets import load_breast_cancer
from causallearn.utils.cit import CIT
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, train_test_split


def stat_test(D, V_idx, T_idx, S=None):
    kci_obj = CIT(D, "kci")  # construct a CIT instance with data and method name
    pValue = kci_obj(V_idx, T_idx, S)
    return pValue


if __name__=='__main__':
    D = load_breast_cancer(as_frame=True).frame.values
    V_indices = np.arange(D.shape[1] - 1)
    T_idx = D.shape[1] - 1


def forward_selection(D, V_indices, T_idx, a):
    # First, we will initialize an empty set for selected variables
    S = set()
    
    # And we will also initialize the set of remaining variables (subtracting the selected variables (S) from all variables (V_indices))
    R = set(V_indices) - S

    # We will iterate until there are no changes in the set of selected variables (S)
    while True:
        # Store a copy of S to check later if it changes
        S_prev = S.copy()

        # Initialize the minimum p-value to infinity and V_star to None
        min_p_value = float('inf')
        V_star = None # Variable currently considered for inclusion

        # We will check each variable in the R set
        for V_idx in R:
            # Calculate the p-value for the current variable against the target variable, conditioning on S
            p_value = stat_test(D, V_idx, T_idx, list(S))

            # If the p-value is lower than the current minimum, update the minimum and remember the variable
            if p_value < min_p_value:
                min_p_value = p_value
                V_star = V_idx

        # Now we can check if a variable was identified
        if V_star is not None:
            # Remove the identified variable from the remaining set (R)
            R.remove(V_star)

            # We check if the variable under consideration for includsion meets our criterion of having a p-val lower than our threshold alpha
            if min_p_value <= a:
                # And if it does, we add it to our selected variables
                S.add(V_star)

        # If S hasn't changed in this iteration, break the loop (while S changes...)
        if S == S_prev:
            break

    return list(S)


def backward_selection(D, T_idx, S, a):
    # Convert S to a set if it's not already one
    S = set(S)

    # Again, we will iterate until there are no changes in the set of selected variables (S)
    while True:
        S_prev = S.copy()  # to check if S changes
        max_p_value = -float('inf')
        V_star = None

        # This time, we iterate over each variable in the selected set (S)
        for V_idx in S:
            S_without_V = S - {V_idx}
            # Calculate the p-value for the current variable against the target variable, conditioning on S without the current variable
            p_value = stat_test(D, V_idx, T_idx, list(S_without_V))

            # If the p-value is higher than the current maximum, update the maximum and remember this variable
            if p_value > max_p_value:
                max_p_value = p_value
                V_star = V_idx

        # Again the idea is that we check if a variable was identified
        if V_star is not None:
            # and then if the p-value of the identified variable is greater than the threshold, we remove it from S 
            if max_p_value > a:
                S.remove(V_star)

        # And this is our condition as previously to break the while loop (while S changes..)
        if S == S_prev:
            break
        
    return list(S)


def cross_validate_model(X, y, classifier, params, k_folds=5):
    # Initialize a StratifiedKFold object for k-fold cross-validation, because stratification ensures each fold is a good representative of the whole
    # We will shuffle data as we did in Assignment 6 before splitting into batches
    # Using random_state here to check for integrity of results while testing the functions. Will reproduce the results, unlike on assignment 6 where random_state = None
    skfold = StratifiedKFold(n_splits=k_folds, shuffle=True, random_state=42)

    # List to store the ROC AUC scores for each fold
    roc_auc_scores = []

    # Iterate over each fold and keep the training and validation data (test)
    for train_index, test_index in skfold.split(X, y):
        X_train, X_val = X[train_index], X[test_index]
        y_train, y_val = y[train_index], y[test_index]

        # Initialize the classifier with the given parameters. We will use the same format as Assignment 6
        model = classifier(**params) # '**params' unpacks the dictionary into the function arguments

        # Train the model
        model.fit(X_train, y_train)

        # Predict probabilities on the validation data
        y_pred = model.predict_proba(X_val)[:, 1] # 'predict_proba' outputs probabilities for each class and we select those for the positive  class

        # ROC AUC score for this fold
        roc_auc = roc_auc_score(y_val, y_pred)
        roc_auc_scores.append(roc_auc)

    # Calculate the mean of the ROC AUC scores from all folds to get an overall performance metric for the model
    return np.mean(roc_auc_scores)


# We will create some plots to compare cross-fold validation
def plot_scores(scores, configs, title):
    best_score_index = scores.index(max(scores))  # Index of the best score
    colors = ['blue' if i != best_score_index else 'red' for i in range(len(scores))]
    plt.bar(configs, scores, color=colors)
    plt.xlabel('Configuration')
    plt.ylabel('ROC AUC Score')
    plt.title(title)
    plt.show()

# and also the roc curve compared to the naive classifier as we did in Assignment 6
def plot_roc_curve(fpr, tpr, roc_auc, model_name):
    plt.figure()
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label = 'Naive Classifier')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curve of {model_name}')
    plt.legend(loc="lower right")
    plt.show()




### Starting the analysis ###
features = D[:, :-1]
target_var = D[:, -1]

# Splitting the data into train (80%) and test sets (20%) and stratifying based on the target variable to maintain the target distribution
X_train, X_test, y_train, y_test = train_test_split(features, target_var, test_size = 0.2, stratify = target_var, random_state = 42)

# We need to merge the X train and Y train, since our forward and backward selection use the dataframe with the included target feature, (T_idx)
D_train = np.concatenate([X_train, y_train.reshape(-1, 1)], axis=1)

# Let us apply forward and backward selection, using the 3 different thresholds of alpha
a = [0.05, 0.01, 0.005]
forward_res = []
backward_res = []
count = 0
for thresh in a:
    for_sel = forward_selection(D = D_train, V_indices = V_indices, T_idx = T_idx, a = thresh)
    back_sel = backward_selection(D = D_train, T_idx = T_idx, S = V_indices, a = thresh)
    forward_res.append(for_sel)
    backward_res.append(back_sel)
    print(f'For a = {thresh}:\nForward_Selection_Variables: {forward_res[count]}\nBackward_Selection_Variables: {backward_res[count]}\n\n')
    count += 1


# Filter X_train and X_test to include only the selected features. All the forward selection results are the same, so we can use any one.
# We could do the same with backward selection but we will 
selected_features = forward_res[0]
X_train_selected = X_train[:, selected_features]
X_test_selected = X_test[:, selected_features]


############### Using Selected Features ###############

# SVM configurations
# L2 penalty is default for 'linear'
#svm_params = [ {'kernel': 'linear', 'C': c} for c in [0.1, 1, 5, 10, 20] ] 
svm_params = [{'kernel': 'linear', 'C': c, 'probability': True} for c in [0.1, 1, 5, 10, 20]]
svm_scores = []
svm_config_labels = []


# Random Forest configurations
rf_params = [ {'n_estimators': n, 'criterion': 'entropy'} for n in [25, 50, 100, 150, 200] ]
rf_scores = []
rf_config_labels = []


# Evaluate SVM configurations
for params in svm_params:
    average_roc_auc = cross_validate_model(X_train_selected, y_train, SVC, params, k_folds=5)
    print("SVM with params {}: Average ROC AUC: {}".format(params, average_roc_auc))
    svm_scores.append(average_roc_auc)
    svm_config_labels.append('C={}'.format(params['C']))


# Evaluate Random Forest configurations
for params in rf_params:
    average_roc_auc = cross_validate_model(X_train_selected, y_train, RandomForestClassifier, params, k_folds=5)
    print("Random Forest with params {}: Average ROC AUC: {}".format(params, average_roc_auc))
    rf_scores.append(average_roc_auc)
    rf_config_labels.append('n={}'.format(params['n_estimators']))

# While there is no significant reason based on the results to plot, let us do it anyway
plot_scores(svm_scores, svm_config_labels, 'SVM ROC AUC Scores (linear - L2)')
plot_scores(rf_scores, rf_config_labels, 'Random Forest ROC AUC Scores (entropy)')


# We will keep the best configurations that we got from the cross fold validation for both SVM and RF
best_svm_params = {'kernel': 'linear', 'C': 20} 
best_rf_params = {'n_estimators': 50, 'criterion': 'entropy'} 

# Train the best SVM model on the full training set
best_svm_model = SVC(**best_svm_params, probability=True)
best_svm_model.fit(X_train_selected, y_train)

# Evaluate the best SVM model on the test set
svm_predictions = best_svm_model.predict_proba(X_test_selected)[:, 1]
svm_roc_auc = roc_auc_score(y_test, svm_predictions)
print("SVM: ROC AUC on test set:", svm_roc_auc)

# Train the best Random Forest model on the full training set
best_rf_model = RandomForestClassifier(**best_rf_params)
best_rf_model.fit(X_train_selected, y_train)

# Evaluate the best Random Forest model on the test set
rf_predictions = best_rf_model.predict_proba(X_test_selected)[:, 1]
rf_roc_auc = roc_auc_score(y_test, rf_predictions)
print("Random Forest: ROC AUC on test set:", rf_roc_auc)

# Calculate ROC Curve for SVM
fpr_svm, tpr_svm, _ = roc_curve(y_test, svm_predictions)
roc_auc_svm = auc(fpr_svm, tpr_svm)

# Calculate ROC Curve for Random Forest
fpr_rf, tpr_rf, _ = roc_curve(y_test, rf_predictions)
roc_auc_rf = auc(fpr_rf, tpr_rf)

plot_roc_curve(fpr_svm, tpr_svm, roc_auc_svm, 'SVM')
plot_roc_curve(fpr_rf, tpr_rf, roc_auc_rf, 'Random Forest')


print(f'Selected features: {forward_res[0]}')
breast_cancer_data = load_breast_cancer()
feature_names = breast_cancer_data.feature_names
print("Feature names:", feature_names[forward_res])


############### Using all of the Features ###############
### To re-run without feature selection, we can just use the original x_train, without filtering the features. Basicly, we will do the same
### analysis by changing the variables from X_train_selected (filtered features) to X_train (original all-features)


svm_scores_all = []
svm_config_labels_all = []

rf_scores_all = []
rf_config_labels_all = []

for params in svm_params:
    average_roc_auc_all = cross_validate_model(X_train, y_train, SVC, params, k_folds=5)
    print("SVM with params {}: Average ROC AUC: {}".format(params, average_roc_auc_all))
    svm_scores_all.append(average_roc_auc_all)
    svm_config_labels_all.append('C={}'.format(params['C']))


for params in rf_params:
    average_roc_auc_all = cross_validate_model(X_train, y_train, RandomForestClassifier, params, k_folds=5)
    print("Random Forest with params {}: Average ROC AUC: {}".format(params, average_roc_auc_all))
    rf_scores_all.append(average_roc_auc_all)
    rf_config_labels_all.append('n={}'.format(params['n_estimators']))

plot_scores(svm_scores_all, svm_config_labels_all, 'SVM ROC AUC Scores, no selection (linear - L2)')
plot_scores(rf_scores_all, rf_config_labels_all, 'Random Forest ROC AUC Scores, no selection (entropy)')


best_svm_params_all = {'kernel': 'linear', 'C': 5} 
best_rf_params_all = {'n_estimators': 100, 'criterion': 'entropy'} 

# Train the best SVM model on the full training set
best_svm_model_all = SVC(**best_svm_params_all, probability=True)
best_svm_model_all.fit(X_train, y_train)

# Evaluate the best SVM model on the test set
svm_predictions_all = best_svm_model_all.predict_proba(X_test)[:, 1]
svm_roc_auc_all = roc_auc_score(y_test, svm_predictions_all)
print("SVM: ROC AUC on test set:", svm_roc_auc_all)

# Train the best Random Forest model on the full training set
best_rf_model_all = RandomForestClassifier(**best_rf_params_all)
best_rf_model_all.fit(X_train, y_train)

# Evaluate the best Random Forest model on the test set
rf_predictions_all = best_rf_model_all.predict_proba(X_test)[:, 1]
rf_roc_auc_all = roc_auc_score(y_test, rf_predictions_all)
print("Random Forest: ROC AUC on test set:", rf_roc_auc_all)

# Calculate ROC Curve for SVM
fpr_svm_all, tpr_svm_all, _ = roc_curve(y_test, svm_predictions_all)
roc_auc_svm_all = auc(fpr_svm_all, tpr_svm_all)

# Calculate ROC Curve for Random Forest
fpr_rf_all, tpr_rf_all, _ = roc_curve(y_test, rf_predictions_all)
roc_auc_rf_all = auc(fpr_rf_all, tpr_rf_all)

plot_roc_curve(fpr_svm_all, tpr_svm_all, roc_auc_svm_all, 'SVM')
plot_roc_curve(fpr_rf_all, tpr_rf_all, roc_auc_rf_all, 'Random Forest')


#### To run forward backward selection on all data, we will use the original dataset -- D ###

for_sel_all = forward_selection(D = D, V_indices = V_indices, T_idx = T_idx, a = 0.05)
back_sel_all = backward_selection(D = D, T_idx = T_idx, S = for_sel_all, a = 0.05)
forward_res.append(for_sel)
backward_res.append(back_sel)
print(f'For a = 0.05:\nForward_Selection_Variables: {for_sel_all}\nBackward_Selection_Variables: {back_sel_all}\n\n')

