import tempfile
import unittest
from pathlib import Path
import numpy as np
from analysis import gaussian, fit_profile, generate_csv, read_csv

class Tests(unittest.TestCase):
    def test_exact_width(self):
        x=np.linspace(-3,3,201)
        m,_=fit_profile(x,gaussian(x,2,.2,.6,.1))
        self.assertAlmostEqual(m['sigma_mm'],.6,places=7)
        self.assertAlmostEqual(m['center_mm'],.2,places=7)
        self.assertAlmostEqual(m['fwhm_mm'],2*np.sqrt(2*np.log(2))*.6,places=7)
    def test_csv_batch(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'sample.csv'; generate_csv(path)
            runs=read_csv(path); self.assertEqual(len(runs),3)
            for (name,data), sigma in zip(runs.items(),[.4,.55,.7]):
                metrics,_=fit_profile(data[:,0],data[:,1])
                self.assertLess(abs(metrics['sigma_mm']-sigma),.02)
    def test_flat(self):
        with self.assertRaises(ValueError): fit_profile(np.arange(20),np.ones(20))
    def test_bad_positions(self):
        with self.assertRaises(ValueError): fit_profile(np.zeros(20),np.arange(20))

if __name__ == '__main__': unittest.main()
