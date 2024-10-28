import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# read data
data = pd.read_csv("egfr.txt",delimiter="\t",header=None,)

# Create a histogram
plt.hist(data.iloc[:, 2] , bins=25, alpha=0.5)
plt.xlim(0, 20)
# Add title and labels
plt.title('Histogram of EGFR loci')
plt.xlabel('Read depth')
plt.ylabel('Frequency')

# Save the plot
plt.savefig('egfr.png')

# Optional: Show the plot
plt.show()

