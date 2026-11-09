import unittest
import numpy as np
from propagation import arrival_times,front_points,invert


class PropagationTests(unittest.TestCase):
    def setUp(self):
        self.stations=np.array([[-30,-25],[0,-30],[30,-25],[-35,0],[0,0],[35,0],[-30,25],[0,30],[30,25]],float)
        self.source=np.array([-12.,-8.])

    def test_circle_without_wind_is_distance_over_speed(self):
        expected=-60+np.linalg.norm(self.stations-self.source,axis=1)/(6*.06)
        np.testing.assert_allclose(arrival_times(self.stations,self.source,-60,6),expected)

    def test_circle_rotation_and_ellipse_ratio_one_are_equivalent(self):
        a=arrival_times(self.stations,self.source,-60,6,(2,1),1,0)
        for angle in [30,75,150]:
            np.testing.assert_allclose(a,arrival_times(self.stations,self.source,-60,6,(2,1),1,angle))

    def test_station_is_on_elliptical_front_at_its_arrival(self):
        speed=6;ratio=1.6;angle=30;wind=np.array([2,1])
        times=arrival_times(self.stations,self.source,-60,speed,wind,ratio,angle)
        theta=np.deg2rad(angle);rotation=np.array([[np.cos(theta),np.sin(theta)],[-np.sin(theta),np.cos(theta)]])
        for station,time in zip(self.stations,times):
            elapsed=time+60
            d=(station-self.source-wind*.06*elapsed)@rotation.T
            value=(d[0]/(speed*.06*np.sqrt(ratio)*elapsed))**2+(d[1]/(speed*.06/np.sqrt(ratio)*elapsed))**2
            self.assertAlmostEqual(value,1.,places=8)

    def test_no_front_before_release_and_linear_growth(self):
        self.assertEqual(len(front_points(self.source,-60,-61,6)),0)
        a=front_points(self.source,-60,-40,6)-self.source
        b=front_points(self.source,-60,-20,6)-self.source
        np.testing.assert_allclose(b,2*a,atol=1e-12)

    def test_inverse_recovers_known_source_on_distributed_network(self):
        observed=arrival_times(self.stations,self.source,-60,6)
        best,_=invert(self.stations,observed,starts=12)
        self.assertLess(np.hypot(best['x']+12,best['y']+8),.01)
        self.assertAlmostEqual(best['speed'],6,places=3)
        self.assertLess(best['rmse_minutes'],1e-5)

    def test_fixed_speed_has_three_unknowns(self):
        observed=arrival_times(self.stations,self.source,-60,6)
        best,_=invert(self.stations,observed,fixed_speed=6,starts=5)
        self.assertEqual(best['parameters'],3)
        self.assertEqual(best['jacobian_rank'],3)

    def test_supercritical_fixed_scenario_is_reported_not_silently_fitted(self):
        best,solutions=invert(self.stations,np.arange(9),wind=(8,0),fixed_speed=3)
        self.assertIsNone(best);self.assertEqual(solutions,[])


if __name__=='__main__':unittest.main()
