#Before using this script, make sure backend is set to "tensorflow" in $HOME/.keras/keras.json and make sure that gpu_keras_theano environment is activated

from tensorflow.keras.models import load_model
import numpy as np

monocyte_regression = load_model("/projects/pfenninggroup/machineLearningForComputationalBiology/eramamur_stuff/ml_monocyte/model_8.hdf5", compile = False)
monocyte_regression.save("/projects/pfenninggroup/machineLearningForComputationalBiology/eramamur_stuff/ml_monocyte/model_8_tf.hdf5")