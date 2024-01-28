import pyutil.util as pu
"""
Check following 3 subroutines in 'ext/feedforward.f90'
 - subroutine ff_change_activation()
 - subroutine ff_print_info()
 - element subroutine ff_activate()

    ActFunc     function type               Indexing
    ------------------------------------------------
    linear      Linear                      0
    tanh        hypoerbolic tan             1
    sigmoid     sigmoid                     2
    mtanh       scaled tanh                 3
    twist       mtanh + linear twisting     4
    relu        ReLU                        5

"""

def activate(input_layer, iActFunc):
    from tensorflow import nn as tfnn
    if iActFunc == 0:
        return input_layer
    elif iActFunc == 1:
        output_layer = tfnn.tanh(input_layer, name='tanh')
    elif iActFunc == 2:
        output_layer = tfnn.sigmoid(input_layer, name='sigmoid')
    elif iActFunc == 5:
        output_layer = tfnn.relu(input_layer, name='relu')
    else:
        pu.call_error('input error - activation.activate',
                'not supported function type - {}'.format(iActFunc))
        exit(1)
    return output_layer

def conv_to_idx(ActFunc):
    """
conv_to_idx(ActFunct)
------------------------------------------------------------------
Convert type of activation function to corresponding integer index.
This function should be compatible with FEEDFORWARD function of 
AENET SOURCE for prediction program.
    """
    if 'str' not in str(type(ActFunc)):
        pu.call_error('input error - activation.indexing_activation_functions',
            'string input only - {}'.format(type(ActFunc)))
    if ActFunc in ('linear','l'):  
        return 0
    elif ActFunc in ('tanh','t'):    
        return 1
    elif ActFunc in ('sigmoid','s'): 
        return 2
    elif ActFunc in ('mtanh','m'):   
        return 3
    elif ActFunc in ('twist','w'):
        return 4
    elif ActFunc in ('relu'):
        return 5
    else:
        pu.call_error('invalid input warning - activation.indexing_activation_functions',
                '\'{}\' is not supported function type - set to linear function'.format(ActFunc),
                isexit=False)
        return 0

def conv_to_str(ActFunc):
    """
conv_to_str(ActFunct)
------------------------------------------------------------------
Convert type of activation function to corresponding name.
This function should be compatible with FEEDFORWARD function of 
AENET SOURCE for prediction program.
    """
    if 'int' not in str(type(ActFunc)):
        pu.call_error('input error - activation.indexing_activation_functions',
            'integer input only - {}'.format(type(ActFunc)))
    if ActFunc == 0:
        return 'linear'
    elif ActFunc == 1:
        return 'tanh'
    elif ActFunc == 2:
        return 'sigmoid'
    elif ActFunc == 3:
        return 'mtanh'
    elif ActFunc == 4:
        return 'twist'
    elif ActFunc == 5:
        return 'relu'
    else:
        pu.call_error('invalid input warning - activation.indexing_activation_functions',
                '\'{}\' is not supported function type - set to linear function'.format(ActFunc),
                isexit=False)
        return 'linear'
