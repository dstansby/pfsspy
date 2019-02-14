# Write the benchmarking functions here.
# See "Writing benchmarks" in the asv docs for more information.
# import pfsspy
import numpy as np
import pfsspy


class TimeSuite:
    """
    An example benchmark that times the performance of various kinds
    of iterating over dictionaries in Python.
    """
    def time_dipole(self):
        ntheta = 30
        nphi = 20
        nr = 10
        rss = 2.5

        phi = np.linspace(0, 2 * np.pi, nphi)
        theta = np.linspace(-np.pi / 2, np.pi / 2, ntheta)
        theta, phi = np.meshgrid(theta, phi)

        def dipole_Br(r, theta):
            return 2 * np.sin(theta) / r**3

        br = dipole_Br(1, theta).T
        input = pfsspy.Input(br, nr, rss)
        output = pfsspy.pfss(input)


'''
class TimeSuite:
    """
    An example benchmark that times the performance of various kinds
    of iterating over dictionaries in Python.
    """
    def setup(self):
        self.d = {}
        for x in range(500):
            self.d[x] = None

    def time_keys(self):
        for key in self.d.keys():
            pass

    def time_iterkeys(self):
        for key in self.d.iterkeys():
            pass

    def time_range(self):
        d = self.d
        for key in range(500):
            x = d[key]

    def time_xrange(self):
        d = self.d
        for key in xrange(500):
            x = d[key]


class MemSuite:
    def mem_list(self):
        return [0] * 256
'''
