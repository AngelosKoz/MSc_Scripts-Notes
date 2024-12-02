# %%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.svm import LinearSVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve, auc
from sklearn.calibration import CalibratedClassifierCV

A_df = (pd.read_csv('Dataset6.A_XY.csv', header = None))
A_X = A_df.values[:, :-1]
A_Y = A_df.values[:, -1].ravel()
A_df_X = A_df.iloc[:, :-1]  # All columns except the last one
A_df_Y = A_df.iloc[:, -1]   # The last column



##### PART 1 - Data Exploration #####

# Counting the number of unique values and analyzing the range (min/max) of these values can give us insight to 
# distinguish between categorical and continuous variables

# 1) Count Unique Values for each feature (column)
# 2) Find Minimum and Maximum Values for each feature (column)
# 3) Check if the range (max - min) is close to the count of unique values. We will use a threshold of 1 difference between unique and min max.
# Here we could also do different sorts of classifications for data types (categorical or continuous). We could check difference between the values for example.
# But since this serves as a starting analysis, we can use this arbitrarily and then check the plots and values.
# 4) Classify as Categorical or Continuous:
#   i) If the range is approximately equal to the count of unique values (difference is +/- 1) or if we only have 2 set of unique values, we assign them as categorical
#   ii) If the range significantly exceeds the count of unique values, the feature is likely continuous. 

def find_cat_cont(features):
    variable_types = {}
    for i in range(features.shape[1]):
        feat_col = features[:, i]
        unique_count = len(np.unique(feat_col))
        min_val = np.min(feat_col) # We dont encounter 0 as minimum values, so we dont have to make any adjustments to max - min
        max_val = np.max(feat_col)
        mean_val = np.mean(feat_col)
        variance_val = np.var(feat_col)
    
        # If we only have 2 unique values (sort of binary, but not necessarily 0-1)
        if unique_count == 2:
            variable_types[i] = { 
                                 'Type' : 'Categorical', 
                                 'Unique_Val': unique_count, 
                                 'Max-Min' : (max_val, min_val),
                                 'Mean': mean_val,
                                 'Variance': variance_val
                                 }
        # Range is approximately equal to the count of unique values
        elif max_val - min_val <= unique_count + 1:
            variable_types[i] = { 
                                 'Type' : 'Categorical', 
                                 'Unique_Val': unique_count, 
                                 'Max-Min' : (max_val, min_val),
                                 'Mean': mean_val,
                                 'Variance': variance_val
                                 }
        # Bigger difference of max-min
        else:
            variable_types[i] = { 
                                 'Type' : 'Continuous', 
                                 'Unique_Val': unique_count, 
                                 'Max-Min' : (max_val, min_val),
                                 'Mean': mean_val,
                                 'Variance': variance_val
                                 }
    return variable_types

A_dtype = find_cat_cont(A_X)

# We can check the statistics using the describe() function as they are in pandas df
statistics = A_df_X.describe()
#print(statistics)

# Since we will be using classification algorithms later on, such as SVM and KNN, we would like features that are informative,
# which means we prefer features with high variance, which we can then normalize. Thus we will first do a sorting based on variance.
# After that, we can sort the correlation with target variable using the variances, and keep those with higher correlation, which
# will help us with SVM and Random Forest (since k-means does not use target variable).

# Calculating variance for each feature
variances = A_df_X.var().sort_values(ascending=False)
#print(variances)

# Sorting the correlation with the target variable based on variance order
correlation = A_df_X.corrwith(A_df_Y).loc[variances.index]
sort_correlation = correlation.sort_values(ascending=False)
#print(correlation)
#print(sort_correlation)

# Looking at both sorted based on variances and sorted on correlation, features 4, 12, 0, 10 stand out in the top correlations sorted 
# based on variances as well as having the highest correlation overall (sort_correlation). 
# We will also pick a feature with the lowest negative correlation, which is feature 5. This is also good because we dont seem to have any other negatively correlated features.
# We can check more detailed statistics of those 5 features in our A_dtype, by using their indeces. 

selected_features = [0, 4, 5, 10, 12]

for i in selected_features:
    print(A_dtype[i])

# We picked 3 categorical, and 2 continuous features so we have some variability.

# Initialize statistics table
statistics_table = pd.DataFrame(columns=['Feature', 'Type', 'Unique Values', 'Class Distribution', 'Mean', 'Variance', 'Correlation'])

# We will add them all in one table, so we will create dummy arrays as NaN to add in our df
plt.figure(figsize=(15, 10))
for i, feature in enumerate(selected_features, 1):
    plt.subplot(2, 3, i)
    feat_i = A_dtype[feature]
    print((feat_i['Max-Min']))
    if feat_i['Type'] == 'Continuous':
        unique_val = feat_i['Unique_Val']
        mean_val = feat_i['Mean']
        variance_val = feat_i['Variance']
        plt.hist(A_df[feature], bins = 20, edgecolor = 'black', alpha=0.7)
        plt.xticks(range(int(feat_i['Max-Min'][0])+1)) # While this is not usually appropriate for continuous, after exploration and for these features, it is appropriate.
        # Create NaN to add in table
        class_distribution = np.nan
        
    else:
        unique_val = feat_i['Unique_Val']
        val, counts = np.unique(A_df[feature], return_counts=True)
        class_distribution = dict(zip(val, counts))
        plt.hist(A_df[feature], bins = 20, edgecolor = 'black', alpha=0.7)
        plt.xticks(range(int(feat_i['Max-Min'][0])+1))
        # Create NaN to add in table
        mean_val = np.nan
        variance_val = np.nan
        

    plt.title(f'Histogram of Feature {feature} ({feat_i["Type"]})')
    plt.xlabel('Unique Values')
    plt.ylabel(f'Frequency (Total:{A_df[0].shape[0]})') # Rows 
        
    # Append values to the table
    statistics_table = statistics_table.append({
        'Feature': feature,
        'Type': feat_i['Type'],
        'Unique Values': unique_val,
        'Class Distribution': class_distribution,
        'Mean': mean_val,
        'Variance': variance_val,
        'Correlation': correlation[feature]
    }, ignore_index=True)

plt.tight_layout()
plt.show()

statistics_table
#statistics_table.to_csv('Feature_table.csv',sep = '\t', index=False) # Used this to save the table of selected features


##### PART 2 - Model Selection #####
# For this analysis, i will be selecting SVM, KNN and RF

# We will store the hyperparameters in lists as dictionaries, with each hyperparameter name as a key.

# 1) SVM (Linear Support Vector Machine)
# A) For L1
svm_params_l1 = [
    {'penalty': 'l1', 'C': c, 'dual': False} for c in [0.1, 1, 10]  # L1 penalty
]
# B) For L2
svm_params_l2 = [
    {'penalty': 'l2', 'C': c, 'dual': True} for c in [0.1, 1, 10]     # L2 penalty
]

# 2) KNN (K-Nearest Neighbors)
knn_params = [
    {'n_neighbors': k, 'weights': w} for k in [3, 5, 7] for w in ['uniform', 'distance']
]

# 3) RF (Random Forest)
rf_params = [
    {'n_estimators': n, 'criterion': c} for n in [50, 100, 150] for c in ['gini', 'entropy']
]

# Combining configurations into a single list for easier iteration later on
configurations = []

for config in svm_params_l1 + svm_params_l2:
    configurations.append({'type': 'svm', **config})
##
for config in knn_params:
    configurations.append({'type': 'knn', **config})
##
for config in rf_params:
    configurations.append({'type': 'rf', **config})

configurations


# Create stratified folds
def create_folds(data, k_folds):
    #Splits the data matrix into k stratified folds with shuffling.
    y = data.iloc[:, -1]  # Target variable is in the last column
    skfold = StratifiedKFold(n_splits = k_folds, shuffle = True, random_state = None)
    folds = [test_index for _, test_index in skfold.split(data, y)]
    return folds # List of arrays: Each array contains the indices of the validation set for each fold.


# Run a model with configurations and return model, performance
def run_model(train_data, test_data, configuration):
    # Splitting the training and test data into features (X) and target (y)
    X_train, y_train = train_data.iloc[:, :-1], train_data.iloc[:, -1]
    X_test, y_test = test_data.iloc[:, :-1], test_data.iloc[:, -1]

    # Since we have stored the type of classifier we will use, we can check accordingly and perform the model appropriately.
    # The hyperparameters of each configuration will indexed and assigned to our models accordingly
    if configuration['type'] == 'svm':
        # 1) SVM (Using LinearSVC model)
        model = LinearSVC(penalty=configuration['penalty'], C=configuration['C'], dual=configuration['dual'])
    elif configuration['type'] == 'knn':
        # 2) KNN (Using KNeighborsClassifier model)
        model = KNeighborsClassifier(n_neighbors=configuration['n_neighbors'], weights=configuration['weights'])
    elif configuration['type'] == 'rf':
        # 3) RandomForest (Using RandomForestClassifier model)
        model = RandomForestClassifier(n_estimators=configuration['n_estimators'], criterion=configuration['criterion'])
    
    
    model.fit(X_train, y_train) # Train model
    predictions = model.predict(X_test) # Use test data to predict
    performance = roc_auc_score(y_test, predictions) # ROC AUC score (using sklearn.metrics)
    
    return model, performance


# Estimate the performance of each configuration using K-Fold cross-validation
def performance_estimation(data, validation_indices, configurations_list):
    performances = []  # store performance of each configuration

    # First we iterate through each configuration in the list
    for configuration in configurations_list:
        fold_performances = []  # store performance for each fold

        # Then we iterate through each set of validation indices for cross-validation
        for val_index in validation_indices:
            # We can find indices for the training set, which is all indices not in the validation set
            train_index = np.setdiff1d(np.arange(len(data)), val_index)

            # Now that we found the indices, we can create our train and test data
            train_data, test_data = data.iloc[train_index], data.iloc[val_index]

            # We will now use the run_model function to calculate the performance
            model, fold_performance = run_model(train_data, test_data, configuration)
            fold_performances.append(fold_performance)

        # Calculate and store the average performance across all folds for the current configuration
        avg_performance = np.mean(fold_performances)
        performances.append(avg_performance)

    return performances


# Function to find the best model
def best_model(performances_list, configurations_list):
    # The best model can be found by the index of the best performance (highest ROC AUC score)
    best_index = np.argmax(performances_list)
    return configurations_list[best_index]


# We will try for k_folds = 5 and k_folds = 10
k_folds = [5,10]
for kf in k_folds:
    folds = create_folds(A_df, k_folds = kf)

    # Estimate performance of each configuration
    performances = performance_estimation(A_df, folds, configurations)
    print(f'Performances for k_folds = {kf}:\n{performances}')

    # Find the best configuration
    best_config = best_model(performances, configurations)
    best_performance = max(performances)
    print(f'Best Configuration: {best_config}')
    print(f'Best Performance (ROC AUC): {best_performance}')

    # We can plot all the scores and configurations
    plt.figure(figsize=(10, 6))
    bars = plt.bar(range(len(performances)), performances, color='lightblue')
    plt.xlabel('Configuration')
    plt.ylabel(f'ROC AUC Score (Best: {round(best_performance,2)})')
    plt.title(f'Performance of Different Configurations (K-Folds={kf})')
    plt.axhline(y=best_performance, color='red', linestyle='--', label=f'Best Performance: {best_performance:.2f}')
    plt.xticks(range(len(performances)), [str(config) for config in configurations], rotation=45, ha='right')
    # To highlight the best configuration, we can do the below:
    best_performance_idx = np.argmax(performances)
    bars[best_performance_idx].set_color('green')
    plt.show()

    # Since by creating a table from the original configurations, we will get only a dictionary, we can re-create it in a way
    # that will allow us to create a more personalized table, in the form of the original Table 1
    configurations_tb = []
    for config in svm_params_l1:
        configurations_tb.append({'Classifier': 'SVM', 'Param 1': 'L1', 'Param 2': 'C=' + str(config['C'])})
    for config in svm_params_l2:
        configurations_tb.append({'Classifier': 'SVM', 'Param 1': 'L2', 'Param 2': 'C=' + str(config['C'])})
    for config in knn_params:
        configurations_tb.append({'Classifier': 'KNN', 'Param 1': 'K=' + str(config['n_neighbors']), 'Param 2': 'Weights=' + config['weights']})
    for config in rf_params:
        configurations_tb.append({'Classifier': 'RF', 'Param 1': 'Trees=' + str(config['n_estimators']), 'Param 2': 'Criterion=' + config['criterion']})

    configuration_table = pd.DataFrame(configurations_tb)
    configuration_table['ROC AUC Score'] = performances
    # Setting row titles (index) as configurations 0 to 17
    configuration_table.index = [f"Configuration {i}" for i in range(len(configurations_tb))]
    configuration_table.to_csv(f'k_{kf}_all_feat_configuration_table.csv', sep='\t', index=True) # Used this to save the table of configurations


# Train the model on the whole dataset
# We will use for K-folds = 5, which is SVM, with the configurations: {'type': 'svm', 'penalty': 'l1', 'C': 10, 'dual': False}

best_config = {'type': 'svm', 'penalty': 'l2', 'C': 10, 'dual': True}

if best_config['type'] == 'svm':
    final_model = LinearSVC(penalty=best_config['penalty'], C=best_config['C'], dual=best_config['dual'])
##
elif best_config['type'] == 'knn':
    final_model = KNeighborsClassifier(n_neighbors=best_config['n_neighbors'], weights=best_config['weights'])
##
elif best_config['type'] == 'rf':
    final_model = RandomForestClassifier(n_estimators=best_config['n_estimators'], criterion=best_config['criterion'])

final_model.fit(A_df_X, A_df_Y)

##### PART 3 - Computing the out-of-sample performance, ROC curve #####

B_df = (pd.read_csv('Dataset6.B_XY.csv', header = None))
B_df_X = B_df.iloc[:, :-1]  # All columns except the last one
B_df_Y = B_df.iloc[:, -1]   # The last column

final_predict = final_model.predict(B_df_X)

performance = roc_auc_score(B_df_Y, final_predict)
print(f'Out-of-sample ROC AUC score: {performance}') # We notice that it is a bit higher than the cross-validation.

# Since we need the probability scores for the ROC curve, we need to extract them. Apparently, LinearSVC does not have that option,
# so this following methods ensures we can get the probabilities.
calibrated_svc = CalibratedClassifierCV(LinearSVC(penalty=best_config['penalty'], C=best_config['C'], dual=best_config['dual']))
calibrated_svc.fit(A_df_X, A_df_Y)

# Get probability scores
final_probabilities = calibrated_svc.predict_proba(B_df_X)[:, 1]

# Since our labels are not 0,1 or -1,1 which is required for the ROC curve, we need to do a transformation.
B_df_Y_bin = B_df_Y.replace({1.0: 0.0, 2.0: 1.0}) # We must be careful to keep the order, or else the ROC curve will be reverse, meaning a very small area.

# Calculate FPR and TPR
fpr, tpr, _ = roc_curve(B_df_Y_bin, final_probabilities)
roc_auc = auc(fpr, tpr)

# Plotting the ROC Curve
plt.figure()
plt.plot(fpr, tpr, color='darkorange', lw=2, label=f'ROC curve (area = {roc_auc:.2f})')
plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--', label='Naive Classifier')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('FPR (False Positive Rate)')
plt.ylabel('TPR (True Positive Rate)')
plt.title('ROC Curve')
plt.legend(loc="upper left")
plt.show()


##### PART 4 - Calculating the 95% Confidence intervals #####


bootstrap_iter = 1000
bootstrap_performances = []

for _ in range(bootstrap_iter):
    # First we will bootstrap sample the indices
    bootstrap_indices = np.random.choice(range(len(B_df_Y)), size=len(B_df_Y), replace=True)
    
    # Now we can use them to get the bootstrap samples
    bootstrap_X = B_df_X.iloc[bootstrap_indices]
    bootstrap_Y = B_df_Y.iloc[bootstrap_indices]
    
    # Calculate performance on bootstrap sample
    bootstrap_predictions = final_model.predict(bootstrap_X)
    bootstrap_performance = roc_auc_score(bootstrap_Y, bootstrap_predictions)
    
    # Store the performance
    bootstrap_performances.append(bootstrap_performance)

# Calculate the 95% confidence intervals
lower_bound = np.percentile(bootstrap_performances, 2.5) #2.5th percentile
upper_bound = np.percentile(bootstrap_performances, 97.5) #97.5th percentile


plt.hist(bootstrap_performances, bins=30, alpha=0.7, color='blue')
plt.axvline(lower_bound, color='green', linestyle='dashed', linewidth=2, label=f'2.5th percentile = {lower_bound:.2f}')
plt.axvline(upper_bound, color='red', linestyle='dashed', linewidth=2, label=f'97.5th percentile = {upper_bound:.2f}')
plt.axvline(performance, color='yellow', linewidth=2, label=f'Original Performance = {performance:.2f}')
plt.xlabel('ROC AUC Score')
plt.ylabel('Frequency')
plt.title('Bootstrap Performance Distribution')
plt.legend()
plt.show()
