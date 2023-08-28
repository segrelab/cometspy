import unittest

from cometspy.utils import pick_random_locations, grow_rocks, _find_unoccupied_adjacent, chemostat
from cometspy.layout import layout
from cometspy.params import params

class TestUtils(unittest.TestCase):

    def test_pick_random_locations(self):
        # Use the pick_random_locations function to generate 3 random locations
        # in a 10x10 grid and check that they are all within the grid
        locs = pick_random_locations(3, (0, 10), (0, 10))
        self.assertEqual(len(locs), 3)
        for loc in locs:
            self.assertGreaterEqual(loc[0], 0)
            self.assertLess(loc[0], 10)
            self.assertGreaterEqual(loc[1], 0)
            self.assertLess(loc[1], 10)
        
    def test_grow_rocks(self):
        # Use the grow_rocks function to generate 5 rocks in a 5x5 grid with
        # mean size 3 and check that they are all within the grid
        locs = grow_rocks(5, (0, 5), (0, 5), 3)
        self.assertEqual(len(locs), 5 * 3)
        for loc in locs:
            self.assertGreaterEqual(loc[0], 0)
            self.assertLessEqual(loc[0], 5)
            self.assertGreaterEqual(loc[1], 0)
            self.assertLessEqual(loc[1], 5)
    
    def test_find_unoccupied_adjacent(self):
        # Use the _find_unoccupied_adjacent function to find all unoccupied
        # adjacent locations to a set of occupied locations in a 5x5 grid
        # Starts with the occupied locations (1, 1), (2, 2), (3, 3)
        # The grid looks like this:
        # ----------------------------------------------
        # | (0, 0) | (0, 1) | (0, 2) | (0, 3) | (0, 4) |
        # |    -   |    -   |    -   |    -   |    -   |
        # ----------------------------------------------
        # | (1, 0) | (1, 1) | (1, 2) | (1, 3) | (1, 4) |
        # |    -   |  OCC.  |    -   |    -   |    -   |
        # ----------------------------------------------
        # | (2, 0) | (2, 1) | (2, 2) | (2, 3) | (2, 4) |
        # |    -   |   -    |  OCC.  |    -   |    -   |
        # ----------------------------------------------
        # | (3, 0) | (3, 1) | (3, 2) | (3, 3) | (3, 4) |
        # |    -   |    -   |    -   |  OCC.  |    -   |
        # ----------------------------------------------
        # | (4, 0) | (4, 1) | (4, 2) | (4, 3) | (4, 4) |
        # |    -   |    -   |    -   |    -   |  OCC.  |
        # ----------------------------------------------
        # So the adjacent unoccupied locations are (0, 1), (1, 2), (2, 1), 
        # (3, 4), (4, 3), (2, 3), (1, 0), (3, 2)
        # ----------------------------------------------
        # | (0, 0) | (0, 1) | (0, 2) | (0, 3) | (0, 4) |
        # |    -   |  ADJ.  |    -   |    -   |    -   |
        # ----------------------------------------------
        # | (1, 0) | (1, 1) | (1, 2) | (1, 3) | (1, 4) |
        # |  ADJ.  |  OCC.  |  ADJ.  |    -   |    -   |
        # ----------------------------------------------
        # | (2, 0) | (2, 1) | (2, 2) | (2, 3) | (2, 4) |
        # |    -   |  ADJ.  |  OCC.  |  ADJ.  |    -   |
        # ----------------------------------------------
        # | (3, 0) | (3, 1) | (3, 2) | (3, 3) | (3, 4) |
        # |    -   |    -   |  ADJ.  |  OCC.  |  ADJ.  |
        # ----------------------------------------------
        # | (4, 0) | (4, 1) | (4, 2) | (4, 3) | (4, 4) |
        # |    -   |    -   |    -   |  ADJ.  |  OCC.  |
        # ----------------------------------------------
        # But the function returns a list of 12 locations since it includes
        # duplicates for each of the 3 loci. So we check that the length of
        # the list is 12 and that the set of locations is 8.
        # And then check the locations are not in the occupied locations
        # and that the expected locations are in the list.
        occupied_locs = [(1, 1), (2, 2), (3, 3)]
        unoccupied_adjacent = _find_unoccupied_adjacent(occupied_locs, (0, 5), (0, 5))
        self.assertEqual(len(unoccupied_adjacent), 12)
        self.assertEqual(len(set(unoccupied_adjacent)), 8)
        for loc in unoccupied_adjacent:
            self.assertTrue((loc[0], loc[1]) not in occupied_locs)
        self.assertEqual(set(unoccupied_adjacent),
                         {(0, 1), (1, 2), (2, 1), (3, 4), (4, 3), (2, 3), (1, 0), (3, 2)})
    
    def test_chemostat(self):
        models = [layout(), layout()]  # Replace with actual models
        reservoir_media = {'glc__D_e': 0.01, 'nh4_e': 1000., 'pi_e': 1000.}
        dilution_rate = 0.1
        layout, parameters = chemostat(models, reservoir_media, dilution_rate)
        self.assertAlmostEqual(parameters.all_params['metaboliteDilutionRate'], dilution_rate)
        self.assertAlmostEqual(parameters.all_params['deathRate'], dilution_rate)
        # Add more assertions related to layout and parameters
        
if __name__ == '__main__':
    unittest.main()
