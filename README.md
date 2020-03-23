# Hand Vein Classification (Gender-Age)
* This project aims to build a model that classifies the subjects into female-male and old-young 
 from the images of the subjects' dorsal hand veins.
## Pipeline
* The project can be summarized into 3 main steps:
    * Dataset: Hand vein images of 200 subjects (.. male & .. female// .. old & .. young)
    * Preprocessing: Some processing (filters, extracting ROI) and augmentation techniques
    *  are applied on the images to enhance their quality.
    * Deeplearning: A fine-tuned model that is used to classify the images that is based on VGG-16.

### Dataset
* The Dataset directory has several sub-directories:
* [Enhanced_Augmented_roi dataset:](https://github.com/AmrMahmoud2/HandVeinClassification/tree/master/Dataset/Enhanced_Augmented_roi%20dataset) that has the enhanced images for age and for gender classification.
* [Original_Dataset:](https://github.com/AmrMahmoud2/HandVeinClassification/tree/master/Dataset/Original_Dataset) the original dataset as raw images without preprocessing/division/excluding_prediction_subjects.
* [Prediction_subjects:](https://github.com/AmrMahmoud2/HandVeinClassification/tree/master/Dataset/Prediction%20subjects) images of 4 subjects that were cut out from the original dataset.
* [Raw data (prediction excluded):](https://github.com/AmrMahmoud2/HandVeinClassification/tree/master/Dataset/Raw%20data%20(prediction%20excluded)) Raw images divided into gender-age with excluding the prediction images.
* #### Note: The number of the dataset and the reserved items for prediction can be identified using the "Age_Gender.xlsx" file.

### Preprocessing
* Firstly, The left-hand images were flipped using the [Flipping_Left](https://github.com/AmrMahmoud2/HandVeinClassification/blob/master/Preprocessing/flipping_left.py) file.
* Secondly, enhancement through several filters using the [Enhance_2020](https://github.com/AmrMahmoud2/HandVeinClassification/blob/master/Preprocessing/Enhance_2020.py) file.
* Thirdly, several augmentation techniques are applied to increase the dataset using the [Augmentation_2020](https://github.com/AmrMahmoud2/HandVeinClassification/blob/master/Preprocessing/Augmentation_2020.py) file.
* Fourthly -and lastly-, extracting the region of interest of the images using the [ROI_2020](https://github.com/AmrMahmoud2/HandVeinClassification/blob/master/Preprocessing/ROI_2020.py
) file.
* In case of working with lbp images, then you'll have to use the [lbp](https://github.com/AmrMahmoud2/HandVeinClassification/blob/master/Preprocessing/lbp.py) file.

### Deeplearning
* Two subdirectories can be found: [Age](https://github.com/AmrMahmoud2/HandVeinClassification/tree/master/Neural_Network_Code/Age) and [Gender](https://github.com/AmrMahmoud2/HandVeinClassification/tree/master/Neural_Network_Code/Gender)
* Both contain subdirectories for their "codes" and "results"
* You first need to train the model on the images you're using (either raw images, enhanced or lbp ones) 
* using the train file.
* Then, you'll need to save the model at the desired epoch 
* (depending on the accuracy and loss graphs obtained from training)
* And lastly, you can use this ".h5" saved model for prediction. 
   During the prediction, a    ".csv"   file is made that comprises the image name, the prediction to it and its true label, a confusion matrix can be found at the end of this    ".csv"   file to better assess the prediction.
