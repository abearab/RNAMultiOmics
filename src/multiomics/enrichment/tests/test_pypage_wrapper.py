"""
Unit tests for the pypage_wrapper module.

Tests PAGE (Pathway Analysis of Gene Expression) functionality
using mutual information analysis with pypage when available.
"""

import unittest
import pandas as pd
import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from pypage_wrapper import (
    run_page,
    create_expression_profile_from_dataframe,
    create_gene_sets_from_dict,
    check_pypage_availability,
    get_pypage_info,
    PYPAGE_AVAILABLE
)


class TestPypageWrapper(unittest.TestCase):
    """Test cases for pypage wrapper functionality."""
    
    def setUp(self):
        """Set up test data."""
        np.random.seed(42)
        
        # Create test expression data
        self.genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5', 'GENE6']
        self.samples = ['Sample1', 'Sample2', 'Sample3', 'Sample4']
        
        # Create expression data with some signal
        self.expression_data = pd.DataFrame(
            np.random.normal(0, 1, (len(self.genes), len(self.samples))),
            index=self.genes,
            columns=self.samples
        )
        
        # Add signal to create informative patterns
        # GENE1 and GENE2 high in Sample1 and Sample2
        self.expression_data.loc['GENE1', ['Sample1', 'Sample2']] += 2
        self.expression_data.loc['GENE2', ['Sample1', 'Sample2']] += 1.5
        
        # GENE3 and GENE4 high in Sample3 and Sample4
        self.expression_data.loc['GENE3', ['Sample3', 'Sample4']] += 2
        self.expression_data.loc['GENE4', ['Sample3', 'Sample4']] += 1.5
        
        # GENE5 and GENE6 have mixed patterns
        self.expression_data.loc['GENE5', ['Sample1', 'Sample3']] += 1
        self.expression_data.loc['GENE6'] += np.random.normal(0, 0.1, len(self.samples))
        
        # Create gene sets
        self.gene_sets = {
            'Pathway_A': ['GENE1', 'GENE2'],  # Should be informative
            'Pathway_B': ['GENE3', 'GENE4'],  # Should be informative
            'Pathway_C': ['GENE1', 'GENE3'],  # Mixed pattern
            'Pathway_D': ['GENE5', 'GENE6'],  # Less informative
            'Single_Gene': ['GENE1'],         # Single gene pathway
            'Large_Pathway': ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5'],  # Large pathway
        }
        
        # Small test data for faster tests
        self.small_expression = self.expression_data.iloc[:4, :3].copy()
        self.small_gene_sets = {
            'Set_A': ['GENE1', 'GENE2'],
            'Set_B': ['GENE3', 'GENE4']
        }
    
    def test_check_pypage_availability(self):
        """Test pypage availability check."""
        availability = check_pypage_availability()
        self.assertIsInstance(availability, bool)
        
        # Should match the module-level variable
        self.assertEqual(availability, PYPAGE_AVAILABLE)
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_get_pypage_info(self):
        """Test getting pypage information."""
        info = get_pypage_info()
        
        self.assertIsInstance(info, dict)
        self.assertIn('version', info)
        self.assertIn('available', info)
        self.assertTrue(info['available'])
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_create_expression_profile_from_dataframe(self):
        """Test creating expression profile from DataFrame."""
        profile = create_expression_profile_from_dataframe(self.expression_data)
        
        # Check that we get a valid pypage ExpressionProfile
        from pypage import ExpressionProfile
        self.assertIsInstance(profile, ExpressionProfile)
        
        # Check dimensions
        self.assertEqual(len(profile.gene_names), len(self.genes))
        self.assertEqual(len(profile.sample_names), len(self.samples))
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_create_gene_sets_from_dict(self):
        """Test creating gene sets from dictionary."""
        gene_sets = create_gene_sets_from_dict(self.gene_sets)
        
        # Check that we get a valid pypage GeneSets object
        from pypage import GeneSets
        self.assertIsInstance(gene_sets, GeneSets)
        
        # Check that gene sets are preserved
        self.assertEqual(len(gene_sets.gene_set_names), len(self.gene_sets))
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_run_page_basic(self):
        """Test basic PAGE analysis."""
        results = run_page(
            expression_data=self.expression_data,
            gene_sets=self.gene_sets,
            n_shuffle=10,  # Small number for testing
            alpha=0.05,
            n_bins=5,
            verbose=False
        )
        
        # Check result structure
        self.assertIsInstance(results, pd.DataFrame)
        
        required_cols = ['Term', 'MI', 'P_value', 'Informative', 'Size']
        for col in required_cols:
            self.assertIn(col, results.columns)
        
        # Should have results for all gene sets
        self.assertEqual(len(results), len(self.gene_sets))
        
        # Check data types
        self.assertTrue(pd.api.types.is_numeric_dtype(results['MI']))
        self.assertTrue(pd.api.types.is_numeric_dtype(results['P_value']))
        self.assertTrue(pd.api.types.is_bool_dtype(results['Informative']))
        
        # Check value ranges
        self.assertTrue(all(results['MI'] >= 0))  # MI should be non-negative
        self.assertTrue(all((results['P_value'] >= 0) & (results['P_value'] <= 1)))
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_run_page_parameters(self):
        """Test PAGE analysis with different parameters."""
        # Test with different number of bins
        results_5_bins = run_page(
            self.small_expression,
            self.small_gene_sets,
            n_shuffle=5,
            n_bins=5
        )
        
        results_3_bins = run_page(
            self.small_expression,
            self.small_gene_sets,
            n_shuffle=5,
            n_bins=3
        )
        
        # Both should work
        self.assertEqual(len(results_5_bins), len(self.small_gene_sets))
        self.assertEqual(len(results_3_bins), len(self.small_gene_sets))
        
        # MI values might be different due to different binning
        # but structure should be the same
        self.assertEqual(list(results_5_bins.columns), list(results_3_bins.columns))
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_run_page_alpha_threshold(self):
        """Test PAGE analysis with different alpha thresholds."""
        # Test with strict alpha
        results_strict = run_page(
            self.expression_data,
            self.gene_sets,
            n_shuffle=10,
            alpha=0.01
        )
        
        # Test with relaxed alpha
        results_relaxed = run_page(
            self.expression_data,
            self.gene_sets,
            n_shuffle=10,
            alpha=0.1
        )
        
        # Relaxed should have same or more informative pathways
        informative_strict = sum(results_strict['Informative'])
        informative_relaxed = sum(results_relaxed['Informative'])
        
        self.assertGreaterEqual(informative_relaxed, informative_strict)
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_run_page_single_gene_pathway(self):
        """Test PAGE analysis with single-gene pathways."""
        single_gene_sets = {
            'Single_A': ['GENE1'],
            'Single_B': ['GENE2'],
            'Single_C': ['GENE3']
        }
        
        results = run_page(
            self.expression_data,
            single_gene_sets,
            n_shuffle=5
        )
        
        # Should handle single-gene pathways
        self.assertEqual(len(results), len(single_gene_sets))
        self.assertTrue(all(results['Size'] == 1))
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_run_page_empty_gene_sets(self):
        """Test PAGE analysis with empty gene sets."""
        # Test with completely empty gene sets
        results_empty = run_page(
            self.expression_data,
            {},
            n_shuffle=5
        )
        
        self.assertEqual(len(results_empty), 0)
        
        # Test with gene sets that have no genes in expression data
        non_matching_sets = {
            'NonExistent': ['FAKE_GENE_1', 'FAKE_GENE_2']
        }
        
        results_non_matching = run_page(
            self.expression_data,
            non_matching_sets,
            n_shuffle=5
        )
        
        # Should have results but with 0 size (filtered out genes)
        if len(results_non_matching) > 0:
            self.assertEqual(results_non_matching.iloc[0]['Size'], 0)
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_run_page_invalid_inputs(self):
        """Test PAGE analysis error handling."""
        # Test with invalid expression data
        with self.assertRaises((ValueError, TypeError)):
            run_page(
                "invalid_expression_data",
                self.gene_sets,
                n_shuffle=5
            )
        
        # Test with invalid gene sets
        with self.assertRaises((ValueError, TypeError)):
            run_page(
                self.expression_data,
                "invalid_gene_sets",
                n_shuffle=5
            )
        
        # Test with invalid parameters
        with self.assertRaises(ValueError):
            run_page(
                self.expression_data,
                self.gene_sets,
                n_shuffle=0  # Invalid number of shuffles
            )
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_run_page_reproducibility(self):
        """Test PAGE analysis reproducibility with fixed seed."""
        # This test checks if results are consistent
        # Note: pypage might not have seed control, so this test might be limited
        
        results1 = run_page(
            self.small_expression,
            self.small_gene_sets,
            n_shuffle=5,
            random_state=42
        )
        
        results2 = run_page(
            self.small_expression,
            self.small_gene_sets,
            n_shuffle=5,
            random_state=42
        )
        
        # MI values should be identical (deterministic)
        # P-values might vary due to permutation randomness
        np.testing.assert_array_almost_equal(
            results1['MI'].values,
            results2['MI'].values,
            decimal=10
        )
    
    def test_pypage_not_available_handling(self):
        """Test graceful handling when pypage is not available."""
        if PYPAGE_AVAILABLE:
            self.skipTest("pypage is available, can't test unavailable case")
        
        # When pypage is not available, functions should raise ImportError or return appropriate messages
        with self.assertRaises(ImportError):
            run_page(
                self.expression_data,
                self.gene_sets,
                n_shuffle=5
            )


class TestPypageIntegration(unittest.TestCase):
    """Integration tests for PAGE analysis."""
    
    def setUp(self):
        """Set up integration test data."""
        np.random.seed(123)
        
        # Create larger, more realistic dataset
        n_genes = 50
        n_samples = 8
        
        self.genes = [f'GENE{i:03d}' for i in range(1, n_genes + 1)]
        self.samples = [f'Condition_A_{i}' for i in range(1, 5)] + [f'Condition_B_{i}' for i in range(1, 5)]
        
        # Create expression data with clear patterns
        self.expression_data = pd.DataFrame(
            np.random.normal(0, 0.5, (n_genes, n_samples)),
            index=self.genes,
            columns=self.samples
        )
        
        # Add strong signal patterns
        condition_a_samples = [s for s in self.samples if 'Condition_A' in s]
        condition_b_samples = [s for s in self.samples if 'Condition_B' in s]
        
        # Pathway 1: Upregulated in Condition A
        pathway1_genes = self.genes[:10]
        self.expression_data.loc[pathway1_genes, condition_a_samples] += 2
        
        # Pathway 2: Upregulated in Condition B
        pathway2_genes = self.genes[10:20]
        self.expression_data.loc[pathway2_genes, condition_b_samples] += 2
        
        # Pathway 3: Mixed pattern (some A, some B)
        pathway3_genes = self.genes[20:30]
        self.expression_data.loc[pathway3_genes[:5], condition_a_samples] += 1.5
        self.expression_data.loc[pathway3_genes[5:], condition_b_samples] += 1.5
        
        # Create pathway gene sets
        self.pathways = {
            'Condition_A_Signature': pathway1_genes,
            'Condition_B_Signature': pathway2_genes,
            'Mixed_Signature': pathway3_genes,
            'Random_Genes': self.genes[30:40],
            'Small_Pathway': self.genes[40:43],
            'Large_Pathway': self.genes[:25]
        }
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_page_identifies_informative_pathways(self):
        """Test that PAGE correctly identifies informative pathways."""
        results = run_page(
            self.expression_data,
            self.pathways,
            n_shuffle=20,
            alpha=0.1,
            n_bins=6
        )
        
        # Check that signature pathways have higher MI than random genes
        condition_a_result = results[results['Term'] == 'Condition_A_Signature'].iloc[0]
        condition_b_result = results[results['Term'] == 'Condition_B_Signature'].iloc[0]
        random_result = results[results['Term'] == 'Random_Genes'].iloc[0]
        
        # Condition-specific signatures should have higher MI
        self.assertGreater(condition_a_result['MI'], random_result['MI'])
        self.assertGreater(condition_b_result['MI'], random_result['MI'])
        
        # At least one of the signature pathways should be informative
        signature_informative = (condition_a_result['Informative'] or 
                                condition_b_result['Informative'])
        self.assertTrue(signature_informative)
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_page_pathway_size_effects(self):
        """Test PAGE analysis with different pathway sizes."""
        results = run_page(
            self.expression_data,
            self.pathways,
            n_shuffle=15,
            alpha=0.05
        )
        
        # Check pathway sizes are correctly reported
        small_pathway = results[results['Term'] == 'Small_Pathway'].iloc[0]
        large_pathway = results[results['Term'] == 'Large_Pathway'].iloc[0]
        
        self.assertEqual(small_pathway['Size'], 3)
        self.assertEqual(large_pathway['Size'], 25)
        
        # Larger pathways might have different MI characteristics
        # but both should be processed correctly
        self.assertGreaterEqual(small_pathway['MI'], 0)
        self.assertGreaterEqual(large_pathway['MI'], 0)
    
    @unittest.skipUnless(PYPAGE_AVAILABLE, "pypage not available")
    def test_page_statistical_significance(self):
        """Test PAGE statistical significance assessment."""
        results = run_page(
            self.expression_data,
            self.pathways,
            n_shuffle=50,  # More shuffles for better p-value estimation
            alpha=0.05
        )
        
        # Check that p-values are properly calculated
        for _, row in results.iterrows():
            # P-values should be in valid range
            self.assertGreaterEqual(row['P_value'], 0)
            self.assertLessEqual(row['P_value'], 1)
            
            # Informative flag should match alpha threshold
            if row['P_value'] <= 0.05:
                self.assertTrue(row['Informative'])
            else:
                self.assertFalse(row['Informative'])


if __name__ == '__main__':
    unittest.main()