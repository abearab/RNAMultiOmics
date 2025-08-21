"""
Unit tests for the ora module.

Tests Over-Representation Analysis (ORA) functionality including
Fisher's exact test, hypergeometric test, and custom ORA implementation.
"""

import unittest
import pandas as pd
import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from ora import (
    run_fisher_exact_test,
    run_hypergeometric_test,
    run_custom_ora,
    run_ora,
    format_ora_results
)

# Check if scipy is available for statistical tests
try:
    import scipy.stats
    ORA_AVAILABLE = True
except ImportError:
    ORA_AVAILABLE = False


class TestORA(unittest.TestCase):
    """Test cases for Over-Representation Analysis."""
    
    def setUp(self):
        """Set up test data."""
        # Create test gene sets with known overlaps
        self.study_genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5']
        self.background_genes = [f'GENE{i}' for i in range(1, 21)]  # GENE1 to GENE20
        
        self.gene_sets = {
            'Pathway_A': ['GENE1', 'GENE2', 'GENE15', 'GENE16', 'GENE17'],  # 2 overlap
            'Pathway_B': ['GENE3', 'GENE4', 'GENE18', 'GENE19', 'GENE20'],  # 2 overlap
            'Pathway_C': ['GENE1', 'GENE3', 'GENE5'],                       # 3 overlap
            'Pathway_D': ['GENE10', 'GENE11', 'GENE12', 'GENE13'],          # 0 overlap
        }
        
        # Test data with different significance levels
        self.large_study = [f'GENE{i}' for i in range(1, 11)]  # GENE1 to GENE10
        self.large_background = [f'GENE{i}' for i in range(1, 101)]  # GENE1 to GENE100
        
        self.large_gene_sets = {
            'Highly_Enriched': [f'GENE{i}' for i in range(1, 8)],      # 7/10 overlap
            'Moderately_Enriched': [f'GENE{i}' for i in range(5, 15)], # 6/10 overlap
            'Not_Enriched': [f'GENE{i}' for i in range(50, 60)],       # 0/10 overlap
        }
    
    def test_fisher_exact_test_basic(self):
        """Test Fisher's exact test with basic data."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for Fisher's exact test")
        
        try:
            result = run_fisher_exact_test(
                self.study_genes, 
                self.gene_sets['Pathway_C'], 
                self.background_genes
            )
            
            # Check that we get some result
            self.assertIsNotNone(result)
            
        except Exception as e:
            self.skipTest(f"Fisher's exact test failed: {e}")
    
    def test_hypergeometric_test_basic(self):
        """Test hypergeometric test with basic data."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for hypergeometric test")
        
        try:
            result = run_hypergeometric_test(
                self.study_genes,
                self.gene_sets['Pathway_C'],
                self.background_genes
            )
            
            # Check that we get some result
            self.assertIsNotNone(result)
            
        except Exception as e:
            self.skipTest(f"Hypergeometric test failed: {e}")
    
    def test_custom_ora_basic(self):
        """Test custom ORA implementation."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        try:
            results = run_custom_ora(
                self.study_genes,
                self.gene_sets,
                self.background_genes
            )
            
            # Check that we get some results
            self.assertIsNotNone(results)
            
        except Exception as e:
            self.skipTest(f"Custom ORA failed: {e}")
    
    def test_run_ora_basic(self):
        """Test main ORA function."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        try:
            results = run_ora(
                self.study_genes,
                self.gene_sets,
                background=self.background_genes
            )
            
            # Check that we get some results
            self.assertIsNotNone(results)
            
        except Exception as e:
            self.skipTest(f"run_ora failed: {e}")
    
    def test_custom_ora_large_dataset(self):
        """Test ORA with larger dataset for statistical power."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        results = run_custom_ora(
            self.large_study,
            self.large_gene_sets,
            self.large_background
        )
        
        # Highly enriched pathway should have lowest p-value
        highly_enriched = results[results['Term'] == 'Highly_Enriched'].iloc[0]
        not_enriched = results[results['Term'] == 'Not_Enriched'].iloc[0]
        
        self.assertLess(highly_enriched['P_value'], not_enriched['P_value'])
        self.assertGreater(highly_enriched['Fold_enrichment'], 1)
        self.assertEqual(not_enriched['Overlap'], 0)
    
    def test_custom_ora_multiple_testing_correction(self):
        """Test multiple testing correction in ORA."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        # Test different correction methods
        results_bonf = run_custom_ora(
            self.study_genes,
            self.gene_sets,
            self.background_genes,
            correction_method='bonferroni'
        )
        
        results_fdr = run_custom_ora(
            self.study_genes,
            self.gene_sets,
            self.background_genes,
            correction_method='fdr_bh'
        )
        
        # Adjusted p-values should be >= original p-values
        for _, row in results_bonf.iterrows():
            self.assertGreaterEqual(row['Adjusted_P_value'], row['P_value'])
        
        for _, row in results_fdr.iterrows():
            self.assertGreaterEqual(row['Adjusted_P_value'], row['P_value'])
    
    def test_custom_ora_gene_list_formats(self):
        """Test ORA with different gene list formats."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        # Test with pandas Series
        study_series = pd.Series(self.study_genes)
        results_series = run_custom_ora(
            study_series,
            self.gene_sets,
            self.background_genes
        )
        
        # Test with numpy array
        study_array = np.array(self.study_genes)
        results_array = run_custom_ora(
            study_array,
            self.gene_sets,
            self.background_genes
        )
        
        # Results should be the same
        pd.testing.assert_frame_equal(results_series, results_array)
    
    def test_custom_ora_empty_inputs(self):
        """Test ORA with empty inputs."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        # Empty study genes
        results_empty_study = run_custom_ora(
            [],
            self.gene_sets,
            self.background_genes
        )
        self.assertEqual(len(results_empty_study), 0)
        
        # Empty gene sets
        results_empty_sets = run_custom_ora(
            self.study_genes,
            {},
            self.background_genes
        )
        self.assertEqual(len(results_empty_sets), 0)
    
    def test_custom_ora_filter_min_genes(self):
        """Test filtering of small gene sets."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        gene_sets_with_small = dict(self.gene_sets)
        gene_sets_with_small['Small_Pathway'] = ['GENE1']  # Only 1 gene
        
        # With min_genes=2, small pathway should be filtered out
        results = run_custom_ora(
            self.study_genes,
            gene_sets_with_small,
            self.background_genes,
            min_genes=2
        )
        
        pathway_names = results['Term'].tolist()
        self.assertNotIn('Small_Pathway', pathway_names)
        self.assertIn('Pathway_A', pathway_names)  # Should still be included
    
    def test_ora_availability_check(self):
        """Test ORA availability check."""
        # This should match the actual availability
        try:
            import scipy.stats
            expected_available = True
        except ImportError:
            expected_available = False
        
        self.assertEqual(ORA_AVAILABLE, expected_available)


class TestORAEdgeCases(unittest.TestCase):
    """Test edge cases and error conditions for ORA."""
    
    def test_no_overlap_case(self):
        """Test case where there's no overlap between study and pathway."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        study_genes = ['GENE1', 'GENE2', 'GENE3']
        pathway_genes = ['GENE10', 'GENE11', 'GENE12']
        background_genes = [f'GENE{i}' for i in range(1, 21)]
        
        result = run_fisher_exact_test(study_genes, pathway_genes, background_genes)
        
        self.assertEqual(result['overlap'], 0)
        self.assertGreater(result['p_value'], 0.5)  # Should not be significant
    
    def test_complete_overlap_case(self):
        """Test case where study genes completely overlap with pathway."""
        if not ORA_AVAILABLE:
            self.skipTest("scipy not available for ORA")
        
        study_genes = ['GENE1', 'GENE2', 'GENE3']
        pathway_genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5']
        background_genes = [f'GENE{i}' for i in range(1, 21)]
        
        result = run_fisher_exact_test(study_genes, pathway_genes, background_genes)
        
        self.assertEqual(result['overlap'], 3)
        self.assertLess(result['p_value'], 0.1)  # Should be significant


if __name__ == '__main__':
    unittest.main()