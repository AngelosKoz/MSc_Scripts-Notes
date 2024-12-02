###### Exercise 3 ######
# %%
import numpy as np
import pandas as pd

# X = IxM matrix
# X_dtype = String beeing either "categorical" or "continuous"
# Y: Ix1 vector of classes
# L: Scalar (For Laplace trick/smoothing)
# D_categorical: 1xM vector with D(m) representing # of possible different values 
# that the categorical m can have (ignore for continuous)


### Training Model ###
### a) ###

def train_NBC(X, X_dtype, Y, L=0, D_categorical=None):
    
    #if L == 0:
        #print('Laplace scalar is set to 0 by default.')
    
    nbc_model = {} #Initialize model as a dictionary
    nbc_model['prior'] = {} #Initialize empty prior likelihoods for classes
    
    # Now, let's get the dimensions a.k.a samples and features of the matrix
    I, M =  X.shape #I is assigned to rows (samples) while M is assigned to columns (features)
    classes = np.unique(Y) #Find the number of unique classes
    #print(f'We have {len(classes)} unique classes')
    
    # Now we need to calculate the prior probabilities for each class
    # (the likelihood of a class without accounting for features)
    for cls in classes:
        # We calculate the probabilities while applying the Laplace trick ( We add L to the numerator and L*#classes to the denominator)
        # We simply sum the number of each class occurence, divided by the total samples, after applying Laplace.
        nbc_model['prior'][cls] = (np.sum(Y == cls) + L) / (I + L*len(classes))
        # With this, we know now prior class probabilities, stored in the key 'prior'

    ### CATEGORICAL ###
    ###
    if X_dtype == "categorical":
        #print('Categorical model selected')
        nbc_model['likelihood'] = {} #Initialize empty likelihoods
        
        # Now we need to calculate likelihoods for each feature in our matrix
        # Taking a feature at a time
        for feat in range(M):
            nbc_model['likelihood'][feat] = {}
            
            # We take into account each possible value of each feature
            for fvalue in range(D_categorical[feat]):
                nbc_model['likelihood'][feat][fvalue]={} 
                
                # And now to calculate the likelihoods for each class
                for cls in classes:
                    # First we calculate the numerator of the likelihood for each feature, given a class
                    # (Simplified, we count occurences of each feature for each specific value for each different class)
                    # At the same time, we apply Laplace trick by adding L to the numerator and L*Possible_feature_values to the denominator
                    # Essentially, we solve this: P( Xm = feat_value | Y = y), for m features, given class y
                    numerator = np.sum((X[:, feat] == fvalue) & (Y==cls)) + L
                    #And we do the same for the denominator
                    denominator = np.sum(Y==cls) + L * D_categorical[feat]
                    nbc_model['likelihood'][feat][fvalue][cls] = numerator/denominator
            
    ## Model description for categorical:
    ## The model returned has 2 basic key "nodes", 1 beeing 'prior', which contains a dictionary with each class
    ## as a key, and the prior probability of said class. The second is 'likelihood'.
    ## Likelihood contains 3 layers of dictionaries, starting from each feature/column (0 to M-1), then possible values for that feature,
    ## and lastly the possible class. A pseudo dictionary look like so: 'likelihood': {features : {feature values : {class}}},
    ## each having multiple values depending on features, each feature's possible values, and the possible classes.

                
    ### CONTINUOUS ###
    ###
    if X_dtype == "continuous":
        #print('Continuous model selected')
        #Here instead of likelihoods, we will calculate the parameters, mean and std.
        nbc_model['mean_std'] = {}
        for feat in range(M):
            nbc_model['mean_std'][feat]= {} #Initialize the dictionary for each feature (we calculate column-wise)
            for cls in classes:
                class_values = X[Y == cls, feat] #Keep only the values of the feature from the given class
                class_feat_mean = np.mean(class_values) #Class-specific feature mean
                class_feat_std = np.std(class_values) #Class-specific feature standard deviation
                nbc_model['mean_std'][feat][cls] = (class_feat_mean, class_feat_std)
              
    ## Model description for continuous:
    ## Again the model is a dictionary, with 2 basic key nodes, 1 beeing the same 'prior' likelihoods.
    ## The second one this time is the parameters, essentially mean and standard deviation which is what we need for continuous data.
    ## This time, we have 2 layers, the first beeing again features from 0 to M-1, the second beeing the possible classes.
    ## A pseudo dictionary look like so: 'mean_std': {features : {class}}
    ## Laplace trick/smoothing does not apply to continuous data
    
    return nbc_model


### Prediction Model ###
### b) ###

def predict_NBC(model, X, X_dtype):
    I, M = X.shape # I = samples (rows), M = features (columns)
    sample_predictions = np.zeros(I) #Empty array for class predictions of samples, with a length equal to the samples
    sample_pred_probs = np.zeros(I)
    
    # For the continuous data, we need to compute the PDF (probability density function).
    # We can create another function to simplify this later on.
    # This computes: P(x | y,m )1/sqrt(2π*σ^2)*e^(-(x-μ)^2 / (2*σ^2))
    def get_PDF(x, mean, std):
        exponent = np.exp(-((x - mean) ** 2) / (2 * std ** 2))
        return (1 / (np.sqrt(2 * np.pi) * std)) * exponent
    
    # Now we can start predicting classes for each sample I
    # First we take each sample
    for smpl in range(I):
        class_posterior = {}
        
        # We can use the prior probabilities stored in our model's dictionary (class, likelihood)
        for cls, prior in model['prior'].items():
            
            if X_dtype == "categorical":
                #print('Categorical model selected')
                # We use log (helps with underflow) to calculate the posterior probabilities
                # We want to calculate P(Y = y | Xi), for y class and i samples
                # This can be written as Π[ P(Xi,m = feat_value | Y=y) ] for m features
                post_prob = np.log(prior) + sum(
                    [np.log(
                        model['likelihood'][feat][X[smpl, feat]][cls]
                        ) for feat in range(M)]
                    )
            
            if X_dtype == "continuous":
                #print('Continuous model selected')
                # We have the mean and standard deviation of each feature m for each class y stored in 'mean_std'
                # Given class y follows normal distribution (Gaussian), the likelihood of observing a continuous value x for a feature
                # P(Xi,m = x | Y=y) = 1/sqrt(2π*σ^2)*e^(-(x-μ)^2 / (2*σ^2)) for i,m,y : samples, features, class
    
                post_prob = np.log(prior) + sum( #We sum the pdfs
                    [np.log(
                        get_PDF(
                            x = X[smpl, feat], 
                            mean = model['mean_std'][feat][cls][0], #Input feature mean
                            std = model['mean_std'][feat][cls][1] #Inpute feature std
                            )) for feat in range(M)]) #Loop through all the features
            
            class_posterior[cls] = post_prob
        
        # Extract the maximum probability along with its apointed class
        # We are essentially finding yi = argmax P(Y = y | Xi) for i samples and y class
        predict_class, class_prob = max(class_posterior.items(), key = lambda x: x[1])
        # Keep the predicted class for current sample
        sample_predictions[smpl] = predict_class
        # We can also store the relative (max) class probability
        sample_pred_probs[smpl] = class_prob

    return sample_predictions, sample_pred_probs
            
    
### c) ###
# An easy way to split the data is to use train_test_split from sklearn
from sklearn.model_selection import train_test_split

def analysis_NBC(X, Y, X_dtype, D_categorical=None, L=0, num_runs=100, size=0.25):
    accuracies = [] # Percent of correctly classified samples
    
    #print(f'Number of iterations for accuracy are set to {num_runs}')
    
    #if size == 0.25:
    #    print(f'The default data splitting is set to 75% training and 25% test sets')
    
    # We will repeat the process of training and predicint while counting for the accuracy
    for _ in range(num_runs):
        # Split the data into (75% training set and 25% test set)
        X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=size)
        
        ### Training ###
        if X_dtype == "categorical":
            #print("Categorical model training")
            model = train_NBC(X = X_train, X_dtype = "categorical", Y = Y_train, L = L, D_categorical=D_categorical)
        
        if X_dtype == "continuous":
            #print("Continuous model training")
            model = train_NBC(X = X_train, X_dtype="continuous", Y = Y_train, L = L)
        
        # Use the prediction model, works for both categorical and continuous
        class_predictions = predict_NBC(model=model, X=X_test, X_dtype=X_dtype)
        
        # Calculate accuracy as a proportion (correct predictions / total predictions)
        correct_predictions = np.sum(class_predictions == Y_test)
        total_predictions = len(Y_test)
        accuracy = correct_predictions/total_predictions
        accuracies.append(accuracy)
    
    return np.mean(accuracies)

### Load the data ###
# Load data for categorical as dataframes
X_cat_df = pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_2/Assignment2_Data/DatasetA_X_categorical.csv', header = None)
Y_cat_df = pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_2/Assignment2_Data/DatasetA_Y.csv', header = None)
D_cat_df = pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_2/Assignment2_Data/DatasetA_D_categorical.csv', header = None)
##
#Load data for continuous as dataframes
X_cont_df = pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_2/Assignment2_Data/DatasetB_X_continuous.csv', header = None)
Y_cont_df = pd.read_csv('/home/aggelos/Desktop/Master/Machine_Learning/Tsamardinos/Assignments/Set_2/Assignment2_Data/DatasetB_Y.csv', header = None)
##
# Convert them to numpy arrays
X_cat = X_cat_df.values
Y_cat = Y_cat_df.values.ravel()
D_cat = D_cat_df.values.ravel()
##
X_cont = X_cont_df.values
Y_cont = Y_cont_df.values.ravel()
## Ravel accounts for an issue that occurs with the array dimensions beeing (N,1) instead of (N,)    

# We can create a list of different Laplace scalars to check their influence.
laplace = [0,1,2,4,10,15,30,50,100,1000]


for L in laplace:
    # For categorical data
    print(f'Using Laplace hyperparameter L = {L}:')
    acc_categorical = analysis_NBC(X = X_cat, X_dtype="categorical", Y = Y_cat, D_categorical=D_cat, L=L)
    print(f"Average accuracy for categorical data: {acc_categorical * 100:.2f}%\n")
    
# For continuous data
acc_continuous = analysis_NBC(X_cont, Y_cont, "continuous")
print(f"Average accuracy for continuous data: {acc_continuous * 100:.2f}%\n")


