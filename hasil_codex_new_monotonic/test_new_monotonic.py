import ast
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
import xarray as xr
from coldpool_flexible import Timing, temperature_criterion

class MonotonicChecks(unittest.TestCase):
    def test_only_temperature_and_manifest_changed(self):
        old=ast.parse(Path(__file__).with_name('reference_old_detector.py').read_text())
        new=ast.parse(Path(__file__).with_name('coldpool_flexible.py').read_text())
        a={n.name:ast.dump(n) for n in old.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        b={n.name:ast.dump(n) for n in new.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        self.assertEqual({k for k in a if a[k]!=b[k]},{'temperature_criterion','run'})

    def test_four_mean_increments_keep_raw_drop_tests(self):
        index=pd.date_range('2025-01-01',periods=180,freq='1min')
        values=30.-np.arange(180)*.12
        values[62]+=.15  # A small minute-scale warming within a broad decline.
        da=xr.DataArray(values,dims=['time'],coords={'time':index.to_numpy()})
        mask=temperature_criterion(da,Timing(1))
        for i in range(10,150):
            means=np.array([values[j-2:j+3].mean() for j in [i-5,i,i+5,i+10,i+15]])
            expected=bool(np.all(np.diff(means)<=0) and values[i+15]-values[i-5]<=-2 and values[i]-values[i-5]<=-.5)
            self.assertEqual(bool(mask.values[i]),expected)
        self.assertTrue(mask.values[60])
        self.assertTrue((np.diff(values[55:76])>0).any())

    def test_missing_mean_cannot_pass(self):
        values=30.-np.arange(100)*.2;values[52]=np.nan
        da=xr.DataArray(values,dims=['time'],coords={'time':pd.date_range('2025-01-01',periods=100,freq='1min')})
        self.assertFalse(temperature_criterion(da,Timing(1)).values[50])

if __name__=='__main__': unittest.main()
