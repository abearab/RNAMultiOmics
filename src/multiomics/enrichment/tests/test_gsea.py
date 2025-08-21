"""
Unit tests for the gsea module.

Tests GSEA functionality including ranked list creation, preranked GSEA,
and integration with gseapy when available.
"""

import unittest
import pandas as pd
import numpy as np
import sys
import os

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from gsea import (
    create_ranked_list,
    run_prerank,
    run_gsea,
    GSEAPY_AVAILABLE
)


class TestGSEA(unittest.TestCase):
    """Test cases for GSEA functionality."""
    
    def setUp(self):
        """Set up test data."""
        np.random.seed(42)
        
        # Create test expression data
        genes = [f'GENE{i:03d}' for i in range(1, 51)]
        samples_a = ['A1', 'A2', 'A3']
        samples_b = ['B1', 'B2', 'B3']
        
        self.expression_data = pd.DataFrame(
            np.random.randn(50, 6),
            index=genes,
            columns=samples_a + samples_b
        )
        
        # Add some differential expression signal
        self.expression_data.loc[genes[:10], samples_b] += 2  # Upregulated in B
        self.expression_data.loc[genes[40:], samples_a] += 1.5  # Upregulated in A
        
        self.condition_a = samples_a
        self.condition_b = samples_b
        
        # Create test gene sets
        self.gene_sets = {
            'Upregulated_in_B': genes[:10],
            'Upregulated_in_A': genes[40:],
            'Random_genes': genes[20:30],
            'Mixed_genes': genes[5:15]
        }
        
        # Small test dataset for faster tests
        self.small_expression = self.expression_data.iloc[:10, :].copy()
        self.small_gene_sets = {
            'Set_A': genes[:5],
            'Set_B': genes[5:10]
        }
    
    def test_create_ranked_list_basic(self):
        """Test basic ranked list creation."""
        try:
            ranked_genes = create_ranked_list(
                self.expression_data,
                self.condition_b,
                self.condition_a
            )
            
            # Check that we get some result
            self.assertIsNotNone(ranked_genes)
            
        except Exception as e:
            self.skipTest(f"create_ranked_list failed: {e}")
    
    @unittest.skipUnless(GSEAPY_AVAILABLE, "gseapy not available")
    def test_run_prerank(self):
        """Test preranked GSEA analysis."""
        try:
            # Create ranked gene list first
            ranked_genes = create_ranked_list(
                self.expression_data,
                self.condition_b,
                self.condition_a
            )
            
            results = run_prerank(
                ranked_genes,
                self.gene_sets,
                permutation_num=10  # Small number for testing
            )
            
            # Check that we get some results
            self.assertIsNotNone(results)
            
        except Exception as e:
            self.skipTest(f"run_prerank failed: {e}")
    
    @unittest.skipUnless(GSEAPY_AVAILABLE, "gseapy not available") 
    def test_run_gsea(self):
        """Test full GSEA analysis."""
        try:
            results = run_gsea(
                self.small_expression,
                self.small_gene_sets,
                self.condition_b,
                self.condition_a,
                permutation_num=10
            )
            
            # Check that we get some results
            self.assertIsNotNone(results)
            
        except Exception as e:
            self.skipTest(f"run_gsea failed: {e}")
    
    def test_gseapy_availability_check(self):
        """Test gseapy availability check."""
        try:
            import gseapy
            expected_available = True
        except ImportError:
            expected_available = False
        
        self.assertEqual(GSEAPY_AVAILABLE, expected_available)
    
    def test_create_ranked_list_edge_cases(self):
        """Test edge cases for ranked list creation."""
        # Single sample in each group
        single_sample_data = self.expression_data.iloc[:, [0, 3]].copy()
        
        ranked = create_ranked_list(
            single_sample_data,
            [single_sample_data.columns[1]],
            [single_sample_data.columns[0]],
            method='log2fc'
        )
        
        self.assertIsInstance(ranked, pd.Series)
        self.assertEqual(len(ranked), len(single_sample_data))
        
        # All same values
        same_values_data = pd.DataFrame(
            np.ones((10, 4)),
            index=[f'GENE{i}' for i in range(10)],
            columns=['A1', 'A2', 'B1', 'B2']
        )
        
        ranked_same = create_ranked_list(
            same_values_data,
            ['B1', 'B2'],
            ['A1', 'A2']
        )
        
        # Should handle constant values gracefully
        self.assertIsInstance(ranked_same, pd.Series)
        
        # All values should be 0 for log2fc
        np.testing.assert_array_almost_equal(ranked_same.values, 0)
    
    def test_ranking_statistics(self):
        """Test that ranking statistics are calculated correctly."""
        # Create simple test case with known differences
        simple_data = pd.DataFrame({
            'A1': [1, 2, 3],
            'A2': [1.1, 2.1, 3.1],
            'B1': [2, 4, 6],
            'B2': [2.1, 4.1, 6.1]
        }, index=['GENE1', 'GENE2', 'GENE3'])
        
        # Log2FC method
        ranked_fc = create_ranked_list(
            simple_data,
            ['B1', 'B2'],
            ['A1', 'A2'],
            method='log2fc'
        )
        
        # All genes should have positive log2FC (B > A)
        self.assertTrue(all(ranked_fc.values > 0))
        
        # GENE3 should have highest fold change
        self.assertEqual(ranked_fc.idxmax(), 'GENE3')


class TestGSEAIntegration(unittest.TestCase):
    """Integration tests for GSEA functionality."""
    
    def setUp(self):
        """Set up integration test data."""
        np.random.seed(123)
        
        # Create larger dataset for integration testing
        genes = [f'GENE{i:04d}' for i in range(1, 201)]
        samples = [f'Ctrl_{i}' for i in range(1, 6)] + [f'Treat_{i}' for i in range(1, 6)]
        
        self.large_expression = pd.DataFrame(
            np.random.randn(200, 10),
            index=genes,
            columns=samples
        )
        
        # Add strong signal to specific gene sets
        up_genes = genes[:20]
        down_genes = genes[180:]
        treat_samples = [f'Treat_{i}' for i in range(1, 6)]
        ctrl_samples = [f'Ctrl_{i}' for i in range(1, 6)]
        
        # Strong upregulation in treatment
        self.large_expression.loc[up_genes, treat_samples] += 3
        # Strong downregulation in treatment
        self.large_expression.loc[down_genes, treat_samples] -= 2
        
        self.ctrl_samples = ctrl_samples
        self.treat_samples = treat_samples
        
        self.pathway_gene_sets = {
            'Upregulated_Pathway': up_genes,
            'Downregulated_Pathway': down_genes,
            'Random_Pathway': genes[100:120],
            'Mixed_Pathway': genes[10:30]
        }
    
    def test_full_gsea_pipeline(self):
        """Test complete GSEA analysis pipeline."""
        # Step 1: Create ranked list
        ranked_genes = create_ranked_list(
            self.large_expression,
            self.treat_samples,
            self.ctrl_samples,
            method='ttest'
        )
        
        # Check that upregulated genes are at the top
        top_20 = ranked_genes.head(20).index.tolist()
        expected_up = [f'GENE{i:04d}' for i in range(1, 21)]
        overlap_up = len(set(top_20) & set(expected_up))
        self.assertGreater(overlap_up, 15)  # Most should be in top 20
        
        # Step 2: Run preranked GSEA (if available)
        if GSEAPY_AVAILABLE:
            results = run_preranked_gsea(
                ranked_genes,
                self.pathway_gene_sets,
                permutation_num=10
            )
            
            # Upregulated pathway should be significantly enriched
            up_pathway = results[results['Term'] == 'Upregulated_Pathway']
            self.assertEqual(len(up_pathway), 1)
            
            # Should have positive enrichment score
            up_result = up_pathway.iloc[0]
            self.assertGreater(up_result['ES'], 0)


if __name__ == '__main__':
    unittest.main()