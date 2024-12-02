# %%
import pandas as pd
import numpy as np
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

# First we will create the train function.
# X: The sample data (samples along rows, features along columns).
# Y: The label vector. Y is the class variable we want to predict.
# n trees: the number of trees to grow (Set 100 as default)
# min samples leaf: minimum number of leaf node observations (Set 5 as default)
# We will choose dictionaries since they are fast to iterate through and are easy to understand and interpret, as seen in previous assignments.


def TrainRF(X, Y, n_trees=100, min_samples_leaf=5):
    # Initialize the dictionary.
    # Here we have 2 keys
    # 1) 'trees' : This will include the decision trees within the Random Forest
    # 2) 'features' : This will include the randomly selected features (columns) to train each tree (hyperparameter)
    model = {'trees': [], 'features': []}

    # Now we need the number of features that will be used at each split. We calculate as the square root of total features (columns)
    n_features = int(np.floor(np.sqrt(X.shape[1])))

    # Next, we will create the trees in the forest. The amount of trees depends on the parameter n_trees.
    # Each following step will be done seperately for each tree.
    for _ in range(n_trees):
        # As a first step, we need to create a bootstrapped dataset, which will itroduce variability among the trees
        ind = np.random.choice(X.shape[0], size=X.shape[0], replace=True) # We use replace = True, since bootstrapping uses replacement, unlike permutation
        X_boot = X[ind] # Create the bootstrapped dataset for X
        Y_boot = Y[ind] # And for Y

        # Randomly select a subset of features for this tree (size is equal to square root of total features, which is common heuristic)
        features = np.random.choice(X.shape[1], size=n_features, replace=False)

        # First we will create the tree classifier using the minimum observations at each leaf node
        tree = DecisionTreeClassifier(min_samples_leaf=min_samples_leaf)
        # And now we can train the decision tree with our bootstrapped dataset.
        tree.fit(X_boot[:, features], Y_boot)
        

        # Now we can store the tree and the features used for said tree.
        model['trees'].append(tree)
        model['features'].append(features)

    # While the model returned seems to have values with the same name under "DecisionTreeClassifier(min_samples_leaf=5)]" they are all different trained models
    return model

# Create the prediction function

def PredictRF(model, X):
    # Store the predictions from each tree in a list
    predict = []
    # Now we can combine the trees with their corresponding features and iterate over them
    for tree, features in zip(model['trees'], model['features']):
        # Now we can use the DecisionTreeClassifier predict function to make predictions for the class
        predict.append(tree.predict(X[:, features]).astype(int)) # Store predictions in the list. We use astype(int) because it returns floats.
        
    
    predict = np.array(predict)
    # Determine the final prediction by majority vote across all trees
    # Now since each tree has made different predictions, we will assign the final class based on the majority rule
    # The code below was modified a bit from : https://www.appsloveworld.com/python/3813/how-can-i-apply-np-apply-along-axis-to-combination-of-two-arrays
    # We will have lists on top of lists containing classes for each sample from each tree and we want to find the one that was assigned by the most trees.
    # Bincount will create a list with length equal to the classes, and assign the number it was found (like we did in decision trees for example 1(0) 2(1) etc), but will look like [1,2]
    # Argmax will return the highest of those occurrences. This function is applied on the columns (axis = 0)
    final_predict = np.apply_along_axis(lambda x: np.bincount(x).argmax(), axis=0, arr=predict)
    
    return final_predict




# %% PART B #####
np_DF = (pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_5/Dataset5_XY.csv')).values
df_Y = np_DF[:, -1].ravel()  # Keep the last column of our dataframe, which is the classes Y (ravel makes sure it is 1D array, as seen in previous assignments)
df_X = np_DF[:, :-1]  # Keep all the rest of the columns as our variables (X)

# We will implement the comparisson for min_samples_leaf = 1 and 10 by adding _1 and _10 respectively to the models and code.
# ( Unfortunately i could not fix in time with a loop due to deadline )

# Splitting the dataset into training (70%) and testing (30%) sets
X_train, X_test, Y_train, Y_test = train_test_split(df_X, df_Y, test_size=0.3)

n_trees = 100
min_samples_leaf_1 = 1
rf_model_1 = TrainRF(X = X_train, Y = Y_train, n_trees=n_trees, min_samples_leaf=min_samples_leaf_1) # Create the random forest model
rf_predict_1 = PredictRF(rf_model_1, X_test) # Make predictions based on the model

min_samples_leaf_10 = 10
rf_model_10 = TrainRF(X = X_train, Y = Y_train, n_trees=n_trees, min_samples_leaf=min_samples_leaf_10) # Create the random forest model
rf_predict_10 = PredictRF(rf_model_10, X_test) # Make predictions based on the model

# We can create a function to compute accuracy as well, since it is #correctly_predicted / #total_samples
def get_accuracy(true_class, predict_class):
    correct_pred = sum(total == predict for total, predict in zip(true_class, predict_class))
    return correct_pred / len(true_class)

# Compute accuracies for each tree
# First we create a list to store the accuracy of each tree in the RF. If we used n_trees = 100, then we expect 100 accuracies
tree_accuracies_1 = []
tree_accuracies_10 = []
# Iterate over each tree and its corresponding set of feature indices in the Random Forest model
for tree, features in zip(rf_model_1['trees'], rf_model_1['features']):
    # Use the current tree to make predictions on the test set
    # The predictions are made using only the subset of features that this tree was trained on
    tree_predictions = tree.predict(X_test[:, features])

    # Calculate the accuracy of the current tree.
    # This is done by comparing the tree's predictions (tree_predictions) with the true labels (Y_test)
    accuracy = get_accuracy(Y_test, tree_predictions)
    tree_accuracies_1.append(accuracy)
# Calculate the accuracy of the combined Random Forest (all trees)
# This is done by comparing the Random Forest predictions (rf_predict) with the actual labels (Y_test)
# We expect the accuracy of the Random Forest to be higher than the average accuracy of individual trees since its combines all the trees
rf_accuracy_1 = get_accuracy(Y_test, rf_predict_1)

for tree, features in zip(rf_model_10['trees'], rf_model_10['features']):
    tree_predictions = tree.predict(X_test[:, features])
    accuracy = get_accuracy(Y_test, tree_predictions)
    tree_accuracies_10.append(accuracy)
rf_accuracy_10 = get_accuracy(Y_test, rf_predict_10)





# %%
# Now we will create the histograms
# For min_samples_leaf = 1
plt.hist(tree_accuracies_1, bins=60, alpha=0.7, label='Individual Trees')
plt.axvline(x=np.mean(tree_accuracies_1), color='blue', linestyle='dashed', linewidth=2, label='Mean Tree Accuracy')
plt.axvline(x=rf_accuracy_1, color='red', linestyle='dashed', linewidth=2, label='Random Forest Accuracy')
plt.xlabel(f'Accuracy (Minimum leaf node obs = {min_samples_leaf_1})')
plt.ylabel(f'Number of Trees (total={n_trees})')
plt.title('Individual Tree Accuracies')
plt.legend()
plt.show()

# For min_samples_leaf = 10
plt.hist(tree_accuracies_10, bins=60, alpha=0.7, label='Individual Trees')
plt.axvline(x=np.mean(tree_accuracies_10), color='blue', linestyle='dashed', linewidth=2, label='Mean Tree Accuracy')
plt.axvline(x=rf_accuracy_10, color='red', linestyle='dashed', linewidth=2, label='Random Forest Accuracy')
plt.xlabel(f'Accuracy (Minimum leaf node obs = {min_samples_leaf_10})')
plt.ylabel(f'Number of Trees (total={n_trees})')
plt.title('Individual Tree Accuracies')
plt.legend()
plt.show()

# %%
