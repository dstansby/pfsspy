"""
Re-projecting from CAR to CEA
-----------------------------
"""

from pfsspy import sample_data
from pfsspy import utils
import matplotlib.pyplot as plt

adapt_maps = utils.load_adapt(sample_data.get_adapt_map())
adapt_map_car = adapt_maps[0]
plt.figure()
adapt_map_car.plot()

adapt_map_cea = utils.car_to_cea(adapt_map_car)
plt.figure()
adapt_map_cea.plot()

plt.show()
