import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_pdf(pdf, axis=None, *args, **kwargs):
    if axis is not None:
        axis.plot(pdf.index, pdf.values, *args, **kwargs)
    else:
        plt.plot(pdf.index, pdf.values, *args, **kwargs)