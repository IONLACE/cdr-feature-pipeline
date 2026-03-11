"""
Complete Size-Independent Interface Complementarity Analysis
Reference: Lawrence & Colman (1993), Venkatraman et al. (2009), CAPRI standards

# Basic usage
python interface_analyzer.py receptor.pdb ligand.pdb

# With chain selection
python interface_analyzer.py receptor.pdb ligand.pdb --chain1 A --chain2 B

# Custom parameters
python interface_analyzer.py receptor.pdb ligand.pdb --cutoff 4.5 --output ./my_results
"""

from html import parser
from pathlib import Path
from paths import PATHS
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree, ConvexHull
from scipy.spatial.distance import cdist
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
import warnings
import os
import sys


@dataclass
class InterfaceAtoms:
    """Container for interface atom data"""
    receptor_coords: np.ndarray    # (N, 3) receptor interface atom coordinates
    ligand_coords: np.ndarray      # (M, 3) ligand interface atom coordinates
    receptor_indices: np.ndarray   # Original indices in full structure
    ligand_indices: np.ndarray     # Original indices in full structure
    receptor_residues: List[str]   # Residue identifiers
    ligand_residues: List[str]     # Residue identifiers

@dataclass
class SurfaceData:
    """Container for surface point data"""
    points: np.ndarray      # (K, 3) surface point coordinates
    normals: np.ndarray     # (K, 3) surface normal vectors
    atom_indices: np.ndarray  # (K,) which atom each point belongs to

class InterfaceComplementarityAnalyzer:
    """
    Complete analyzer for size-independent interface complementarity metrics.
    Implements established methods from literature.
    """
    
    def __init__(self, probe_radius: float = 1.4):
        """
        Initialize analyzer.
        
        Parameters:
        -----------
        probe_radius : float
            Solvent probe radius in Å (default 1.4 for water)
        """
        self.probe_radius = probe_radius
        # logger.info(f"Initialized analyzer with probe radius: {probe_radius} Å")
 
    def load_pdb_structure(self, pdb_file: str, 
                          chain_id: Optional[str] = None) -> Dict:
        """
        Load PDB structure and extract atom information.
        
        Returns:
        --------
        dict with 'coords', 'elements', 'residues', 'chain'
        """
        # logger.info(f"Loading PDB structure: {os.path.basename(pdb_file)}")
        
        coords = []
        elements = []
        residues = []
        chains = []
        atom_names = []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    # Check chain
                    current_chain = line[21] if len(line) > 21 else ' '
                    if chain_id and current_chain != chain_id:
                        continue
                    
                    try:
                        # Parse coordinates
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        
                        # Atom name and element
                        atom_name = line[12:16].strip()
                        element = self._infer_element(atom_name)
                        
                        # Residue information
                        res_name = line[17:20].strip()
                        res_num = int(line[22:26])
                        residue_id = f"{res_name}{res_num}"
                        
                        coords.append([x, y, z])
                        elements.append(element)
                        residues.append(residue_id)
                        chains.append(current_chain)
                        atom_names.append(atom_name)
                        
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing line: {line[:60]}... - {e}")
                        continue
        
        if not coords:
            raise ValueError(f"No atoms found in {pdb_file}")
        
        return {
            'coords': np.array(coords),
            'elements': np.array(elements),
            'residues': np.array(residues),
            'chains': np.array(chains),
            'atom_names': np.array(atom_names)
        }
    
    def _infer_element(self, atom_name: str) -> str:
        """Infer element from atom name"""
        # Remove digits
        element = ''.join([c for c in atom_name if not c.isdigit()])
        
        # Common special cases
        special_cases = {
            'CL': 'Cl', 'BR': 'Br', 'NA': 'Na', 'MG': 'Mg',
            'CA': 'Ca', 'FE': 'Fe', 'ZN': 'Zn', 'MN': 'Mn',
            'CU': 'Cu', 'NI': 'Ni'
        }
        
        for key, value in special_cases.items():
            if element.upper().startswith(key):
                return value
        
        # First character is usually the element
        if element:
            return element[0].upper()
        
        return 'C'  # Default
    
    def define_interface(self, receptor_data: Dict, ligand_data: Dict,
                        distance_cutoff: float = 5.0) -> InterfaceAtoms:
        """
        Define interface atoms based on distance cutoff.
        
        Parameters:
        -----------
        receptor_data, ligand_data : dict
            Atom data from load_pdb_structure
        distance_cutoff : float
            Maximum distance to consider atoms as interface (Å)
            
        Returns:
        --------
        InterfaceAtoms object
        """
        # logger.info(f"Defining interface with cutoff: {distance_cutoff} Å")
        
        rec_coords = receptor_data['coords']
        lig_coords = ligand_data['coords']
        
        # Build KD-trees for fast distance calculation
        rec_tree = cKDTree(rec_coords)
        lig_tree = cKDTree(lig_coords)
        
        # Find receptor atoms near ligand
        rec_distances, rec_neighbors = lig_tree.query(rec_coords, 
                                                     distance_upper_bound=distance_cutoff)
        rec_interface_mask = rec_distances < distance_cutoff
        
        # Find ligand atoms near receptor
        lig_distances, lig_neighbors = rec_tree.query(lig_coords,
                                                     distance_upper_bound=distance_cutoff)
        lig_interface_mask = lig_distances < distance_cutoff
        
        # Extract interface atoms
        rec_if_coords = rec_coords[rec_interface_mask]
        lig_if_coords = lig_coords[lig_interface_mask]
        
        # Get original indices
        rec_if_indices = np.where(rec_interface_mask)[0]
        lig_if_indices = np.where(lig_interface_mask)[0]
        
        # Get residue information
        rec_if_residues = receptor_data['residues'][rec_interface_mask].tolist()
        lig_if_residues = ligand_data['residues'][lig_interface_mask].tolist()
        
        return InterfaceAtoms(
            receptor_coords=rec_if_coords,
            ligand_coords=lig_if_coords,
            receptor_indices=rec_if_indices,
            ligand_indices=lig_if_indices,
            receptor_residues=rec_if_residues,
            ligand_residues=lig_if_residues
        )
    
    def generate_interface_surface(self, interface_atoms: InterfaceAtoms,
                                  points_per_atom: int = 20) -> Tuple[SurfaceData, SurfaceData]:
        """
        Generate surface points for interface atoms only.
        
        Parameters:
        -----------
        interface_atoms : InterfaceAtoms
            Interface atom data
        points_per_atom : int
            Number of surface points to generate per atom
            
        Returns:
        --------
        Tuple of (receptor_surface, ligand_surface)
        """

        
        # Generate surface for receptor interface atoms
        rec_surface = self._generate_atom_surface(
            interface_atoms.receptor_coords, points_per_atom
        )
        
        # Generate surface for ligand interface atoms
        lig_surface = self._generate_atom_surface(
            interface_atoms.ligand_coords, points_per_atom
        )
        
        return rec_surface, lig_surface
    
    def _generate_atom_surface(self, atom_coords: np.ndarray,
                             points_per_atom: int) -> SurfaceData:
        """Generate surface points around atoms"""
        if len(atom_coords) == 0:
            return SurfaceData(
                points=np.zeros((0, 3)),
                normals=np.zeros((0, 3)),
                atom_indices=np.array([])
            )
        
        all_points = []
        all_normals = []
        all_atom_indices = []
        
        # Atomic radii (approximate)
        atom_radius = 1.7  # Default carbon radius
        
        # Generate points for each atom
        for i, center in enumerate(atom_coords):
            # Generate points on a sphere
            sphere_points = self._fibonacci_sphere(center, atom_radius, points_per_atom)
            
            # Calculate normals (point outward from atom center)
            normals = (sphere_points - center) / atom_radius
            
            all_points.extend(sphere_points)
            all_normals.extend(normals)
            all_atom_indices.extend([i] * points_per_atom)
        
        points_array = np.array(all_points)
        normals_array = np.array(all_normals)
        indices_array = np.array(all_atom_indices)
        
        # Remove points that are inside other atoms (simplified)
        if len(points_array) > 0:
            tree = cKDTree(atom_coords)
            distances, _ = tree.query(points_array, k=1)
            # Keep points that are close to surface
            surface_mask = np.abs(distances - atom_radius) < 0.5
            points_array = points_array[surface_mask]
            normals_array = normals_array[surface_mask]
            indices_array = indices_array[surface_mask]
        
        return SurfaceData(
            points=points_array,
            normals=normals_array,
            atom_indices=indices_array
        )
    
    def _fibonacci_sphere(self, center: np.ndarray, radius: float,
                         n_points: int) -> np.ndarray:
        """Generate points on a sphere using Fibonacci spiral"""
        indices = np.arange(n_points, dtype=float) + 0.5
        phi = np.arccos(1 - 2 * indices / n_points)
        theta = np.pi * (1 + 5**0.5) * indices
        
        x = np.cos(theta) * np.sin(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(phi)
        
        points = np.stack([x, y, z], axis=1)
        points = points * radius + center
        
        return points
    
    def calculate_lawrence_colman_sc(self, rec_surface: SurfaceData,
                                    lig_surface: SurfaceData,
                                    cutoff: float = 4.0) -> Dict:
        """
        Calculate Lawrence & Colman (1993) Sc statistic for interface.
        
        Parameters:
        -----------
        rec_surface, lig_surface : SurfaceData
            Surface points for receptor and ligand interfaces
        cutoff : float
            Distance cutoff for correlating surface points (Å)
            
        Returns:
        --------
        dict with Sc score and statistics
        """
        
        if len(rec_surface.points) == 0 or len(lig_surface.points) == 0:
            print("No surface points for Sc calculation")
            return {'sc_score': 0.0, 'n_pairs': 0}
        
        # Build KD-trees
        rec_tree = cKDTree(rec_surface.points)
        lig_tree = cKDTree(lig_surface.points)
        
        # Find correlated point pairs
        sc_values = []
        distances = []
        
        # Query receptor points against ligand surface
        rec_dists, lig_indices = lig_tree.query(rec_surface.points, 
                                               distance_upper_bound=cutoff)
        
        for i, (dist, lig_idx) in enumerate(zip(rec_dists, lig_indices)):
            if dist < cutoff:
                dot = np.dot(rec_surface.normals[i], lig_surface.normals[lig_idx])
                # Apply distance weighting (Lawrence & Colman Eq. 3)
                weight = np.exp(-0.5 * (dist**2))
                sc_values.append(dot * weight)
                distances.append(dist)
        
        # Query ligand points against receptor surface (for symmetry)
        lig_dists, rec_indices = rec_tree.query(lig_surface.points,
                                               distance_upper_bound=cutoff)
        
        for j, (dist, rec_idx) in enumerate(zip(lig_dists, rec_indices)):
            if dist < cutoff and dist > 0:  # Avoid double counting
                dot = np.dot(lig_surface.normals[j], rec_surface.normals[rec_idx])
                weight = np.exp(-0.5 * (dist**2))
                sc_values.append(dot * weight)
                distances.append(dist)
        
        if not sc_values:
            print("No correlated surface point pairs found")
            return {'sc_score': 0.0, 'n_pairs': 0}
        
        # Calculate Sc statistic
        sc_raw = np.mean(sc_values)
        # Normalize from [-1, 1] to [0, 1]
        sc_normalized = (sc_raw + 1) / 2
        sc_normalized = max(0.0, min(1.0, sc_normalized))
        
        return {
            'sc_score': float(sc_normalized),
            'sc_raw': float(sc_raw),
            'n_pairs': len(sc_values),
            'mean_distance': float(np.mean(distances)) if distances else 0.0,
            'std_distance': float(np.std(distances)) if len(distances) > 1 else 0.0,
            # 'method': 'Lawrence & Colman (1993)'
        }
    
    def calculate_packing_statistics(self, interface_atoms: InterfaceAtoms,
                                    distance_cutoff: float = 5.0) -> Dict:
        """
        Calculate interface packing statistics.
        
        Parameters:
        -----------
        interface_atoms : InterfaceAtoms
            Interface atom data
        distance_cutoff : float
            Maximum distance for contact consideration (Å)
            
        Returns:
        --------
        dict with packing statistics
        """
        
        rec_coords = interface_atoms.receptor_coords
        lig_coords = interface_atoms.ligand_coords
        
        if len(rec_coords) == 0 or len(lig_coords) == 0:
            return {'error': 'No interface atoms'}
        
        # Calculate all pairwise distances
        all_distances = cdist(rec_coords, lig_coords)
        
        # Find contacts within cutoff
        contact_mask = all_distances < distance_cutoff
        contact_distances = all_distances[contact_mask]
        
        if len(contact_distances) == 0:
            print("No atomic contacts found within cutoff")
            return {
                'n_contacts': 0,
                'contact_density': 0.0,
                'mean_distance': 0.0,
                'min_distance': 0.0
            }
        
        # Calculate statistics
        n_contacts = len(contact_distances)
        n_atoms = len(rec_coords) + len(lig_coords)
        
        # Contact density (contacts per atom)
        contact_density = n_contacts / n_atoms
        
        # Distance statistics
        mean_dist = np.mean(contact_distances)
        std_dist = np.std(contact_distances)
        min_dist = np.min(contact_distances)
        max_dist = np.max(contact_distances)
        
        # Coefficient of variation (uniformity measure)
        cv = std_dist / mean_dist if mean_dist > 0 else 0.0
        
        # Gap volume estimate (simplified)
        # Assuming optimal distance is 3.0Å
        optimal_dist = 3.0
        gap_volume = np.sum((contact_distances - optimal_dist) ** 2) / len(contact_distances)
        
        return {
            'n_contacts': int(n_contacts),
            'contact_density': float(contact_density),
            'mean_distance': float(mean_dist),
            'std_distance': float(std_dist),
            'min_distance': float(min_dist),
            'max_distance': float(max_dist),
            'distance_cv': float(cv),
            'gap_volume_estimate': float(gap_volume),
            'packing_uniformity': float(1.0 / (1.0 + cv)) if cv > 0 else 1.0
        }
    
    def calculate_interface_compactness(self, interface_atoms: InterfaceAtoms) -> Dict:
        """
        Calculate interface compactness measures.
        
        Parameters:
        -----------
        interface_atoms : InterfaceAtoms
            Interface atom data
            
        Returns:
        --------
        dict with compactness statistics
        """
        
        # Combine interface coordinates
        all_if_coords = np.vstack([
            interface_atoms.receptor_coords,
            interface_atoms.ligand_coords
        ])
        
        if len(all_if_coords) < 4:
            return {'convex_hull_volume': 0.0, 'compactness': 0.0}
        
        try:
            # Calculate convex hull
            hull = ConvexHull(all_if_coords)
            hull_volume = hull.volume
            
            # Calculate bounding sphere
            center = all_if_coords.mean(axis=0)
            max_distance = np.max(np.linalg.norm(all_if_coords - center, axis=1))
            sphere_volume = (4/3) * np.pi * (max_distance ** 3)
            
            # Compactness ratio (higher = more compact)
            if sphere_volume > 0:
                compactness = hull_volume / sphere_volume
            else:
                compactness = 0.0
            
            # Surface area to volume ratio (lower = more compact)
            if hull_volume > 0:
                sa_vol_ratio = hull.area / hull_volume
            else:
                sa_vol_ratio = float('inf')
            
            return {
                'convex_hull_volume': float(hull_volume),
                'bounding_sphere_volume': float(sphere_volume),
                'compactness_ratio': float(compactness),
                'surface_area': float(hull.area),
                'surface_to_volume_ratio': float(sa_vol_ratio),
                'n_vertices': len(hull.vertices)
            }
            
        except Exception as e:
            print(f"Convex hull calculation failed: {e}")
            return {'convex_hull_volume': 0.0, 'compactness': 0.0}
    
    def analyze_chemical_complementarity(self, receptor_data: Dict,
                                        ligand_data: Dict,
                                        interface_atoms: InterfaceAtoms) -> Dict:
        """
        Analyze chemical complementarity at interface.
        Simplified version based on residue properties.
        """
        # Residue properties (simplified)
        hydrophobic = {'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'CYS'}
        hydrophilic = {'ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'ASN', 'GLN', 'SER', 'THR', 'TYR'}
        
        # Get unique interface residues
        rec_residues = set(interface_atoms.receptor_residues)
        lig_residues = set(interface_atoms.ligand_residues)
        
        # Classify residues
        rec_hydrophobic = sum(1 for r in rec_residues if r[:3] in hydrophobic)
        rec_hydrophilic = sum(1 for r in rec_residues if r[:3] in hydrophilic)
        
        lig_hydrophobic = sum(1 for r in lig_residues if r[:3] in hydrophobic)
        lig_hydrophilic = sum(1 for r in lig_residues if r[:3] in hydrophilic)
        
        # Calculate complementarity scores
        total_pairs = len(rec_residues) * len(lig_residues) if rec_residues and lig_residues else 1
        
        # Hydrophobic complementarity (hydrophobic-hydrophobic is good)
        hydrophobic_score = (rec_hydrophobic * lig_hydrophobic) / total_pairs
        
        # Charge complementarity (opposite charges are good)
        # Simplified: assume equal positive/negative for now
        charge_score = 0.5  # Placeholder
        
        # Overall chemical complementarity
        chem_score = (hydrophobic_score + charge_score) / 2
        
        return {
            'hydrophobic_score': float(hydrophobic_score),
            'charge_complementarity': float(charge_score),
            'chemical_score': float(chem_score),
            'n_hydrophobic_rec': rec_hydrophobic,
            'n_hydrophilic_rec': rec_hydrophilic,
            'n_hydrophobic_lig': lig_hydrophobic,
            'n_hydrophilic_lig': lig_hydrophilic,
            'total_interface_residues': len(rec_residues) + len(lig_residues)
        }
    
    def run_complete_analysis(self, complex_pdb: str,
                      chain_rec: str, chain_lig: str,
                      distance_cutoff: float = 5.0,
                      output_dir: str = './output') -> Dict:
        """
        Robust single PDB analysis with proper key naming.
        """
        # 1a. Load complex
        print(f"\n1. Loading complex: {os.path.basename(complex_pdb)}")
        complex_data = self.load_pdb_structure(complex_pdb)
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # 1b. Extract chains
        def extract_chain(data, chain_id):
            chain_set = set(str(chain_id).strip())
            if not chain_set:
                raise ValueError(f"No valid chain specification provided: {chain_id}")
            mask = np.isin(data['chains'], list(chain_set))
            if not np.any(mask):
                raise ValueError(f"No atoms found for chain specification: {chain_id}")
            return {
                'coords': data['coords'][mask],
                'elements': data['elements'][mask],
                'residues': data['residues'][mask],
                'chains': data['chains'][mask],
                'atom_names': data['atom_names'][mask]
            }
        
        receptor_data = extract_chain(complex_data, chain_rec)
        ligand_data = extract_chain(complex_data, chain_lig)
        
        # 2. Define interface
        print("2. Defining interface...")
        interface_atoms = self.define_interface(
            receptor_data, ligand_data, distance_cutoff
        )
        
        # 3. Generate interface surfaces
        print("3. Generating interface surfaces...")
        rec_surface, lig_surface = self.generate_interface_surface(interface_atoms)
        
        # 4. Calculate all metrics
        print("4. Calculating complementarity metrics...")
        
        # Lawrence & Colman Sc
        sc_results = self.calculate_lawrence_colman_sc(rec_surface, lig_surface)
        
        # Packing statistics
        packing_stats = self.calculate_packing_statistics(interface_atoms)
        
        # Compactness
        compactness_stats = self.calculate_interface_compactness(interface_atoms)
        
        # Chemical complementarity (simplified)
        chem_stats = self.analyze_chemical_complementarity(
            receptor_data, ligand_data, interface_atoms
        )
        
        # 5. Compile results
        results = {
            'interface_definition': {
                # 'distance_cutoff': distance_cutoff,
                # 'n_receptor_atoms': len(interface_atoms.receptor_coords),
                # 'n_ligand_atoms': len(interface_atoms.ligand_coords),
                'n_receptor_residues': len(set(interface_atoms.receptor_residues)),
                'n_ligand_residues': len(set(interface_atoms.ligand_residues)),
                # 'receptor_residues': list(set(interface_atoms.receptor_residues))[:10],  # First 10
                # 'ligand_residues': list(set(interface_atoms.ligand_residues))[:10]
            },
            'shape_complementarity': sc_results,
            'packing_statistics': packing_stats,
            'compactness': compactness_stats,
            'chemical_complementarity': chem_stats,
        }
        print(f"\n✓ Analysis complete!")
        
        return results
    
    def _create_summary(self, sc_results, packing_stats, 
                       compactness_stats, chem_stats) -> Dict:
        """Create summary statistics"""
        
        # Overall interface quality score (composite)
        weights = {
            'sc': 0.35,
            'packing': 0.30,
            'compactness': 0.20,
            'chemical': 0.15
        }
        
        # Normalize components to [0, 1]
        sc_norm = sc_results['sc_score']
        packing_norm = packing_stats.get('packing_uniformity', 0.5)
        compactness_norm = compactness_stats.get('compactness_ratio', 0.5)
        chem_norm = chem_stats.get('chemical_score', 0.5)
        
        # Calculate weighted score
        overall_score = (
            weights['sc'] * sc_norm +
            weights['packing'] * packing_norm +
            weights['compactness'] * compactness_norm +
            weights['chemical'] * chem_norm
        )
        
        # Quality assessment
        if overall_score > 0.7:
            quality = "Excellent"
        elif overall_score > 0.5:
            quality = "Good"
        elif overall_score > 0.3:
            quality = "Moderate"
        else:
            quality = "Poor"
        
        return {
            'overall_interface_score': float(overall_score),
            'quality_assessment': quality,
            'component_scores': {
                'shape_complementarity': float(sc_norm),
                'packing_uniformity': float(packing_norm),
                'interface_compactness': float(compactness_norm),
                'chemical_complementarity': float(chem_norm)
            }
        }
    
def main():
    """Main function for command-line usage"""
    
    if len(sys.argv) != 2:
        print("Usage: python -m src.analyze_pocket_conservation <pdb_file.pdb>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    pdb_id = Path(pdb_file).stem.split('_')[0]  # Extract PDB ID from filename

    if not Path(pdb_file).is_absolute():
        pdb_file = PATHS.inp_pdb_dir / f"{pdb_id}" / pdb_file

    print(f"Input PDB file: {pdb_file}")

    # parser.add_argument('--cutoff', type=float, default=5.0)
    # parser.add_argument('--probe', type=float, default=1.4)
    probe_radius = 1.4
    cutoff = 5.0

    chain1 = Path(pdb_file).stem.split('_')[2] # ag_chain    
    chain2 = Path(pdb_file).stem.split('_')[1] # nb_chain

    if not chain1 or not chain2:
        print(
            "Could not resolve both chains. Provide --chain1/--chain2, "
            "or use a filename like <pdb>_<chain1>_<chain2>.pdb"
        )
        sys.exit(1)
    

    output_dir = PATHS.meta_dir / Path(pdb_file).stem  / "interface" 

    try:
        # Initialize analyzer
        analyzer = InterfaceComplementarityAnalyzer(probe_radius=probe_radius)
              
        results = analyzer.run_complete_analysis(
            complex_pdb=pdb_file,
            chain_rec=chain1,
            chain_lig=chain2,
            distance_cutoff = cutoff,
            output_dir=str(output_dir)
        )
        if 'shape_complementarity' in results:
            sc = results['shape_complementarity']
            if 'sc_score' in sc:
                print(f"Lawrence & Colman Sc: {sc['sc_score']:.4f}")
            else:
                print(f"Sc results: {sc}")
        else:
            print("Warning: 'shape_complementarity' key not found")
            print(f"Available keys: {list(results.keys())}")
        
        # if 'summary' in results:
        #     summary = results['summary']
        #     if 'overall_interface_score' in summary:
        #         print(f"Overall score: {summary['overall_interface_score']:.4f}")
        #     if 'quality_assessment' in summary:
        #         print(f"Quality: {summary['quality_assessment']}")
        
        # Save results
        os.makedirs(output_dir, exist_ok=True)

        def flatten_for_csv(obj, parent_key='', out=None):
            if out is None:
                out = {}

            if isinstance(obj, dict):
                for key, value in obj.items():
                    new_key = f"{parent_key}_{key}" if parent_key else str(key)
                    flatten_for_csv(value, new_key, out)
            elif isinstance(obj, list):
                if all(not isinstance(item, (dict, list)) for item in obj):
                    out[parent_key] = '|'.join(str(item) for item in obj)
                else:
                    for idx, item in enumerate(obj):
                        new_key = f"{parent_key}_{idx}" if parent_key else str(idx)
                        flatten_for_csv(item, new_key, out)
            else:
                out[parent_key] = obj

            return out

        # Convert numpy types to Python types for flattening
        def convert(obj):
            if isinstance(obj, np.integer):
                return int(obj)
            elif isinstance(obj, np.floating):
                return float(obj)
            elif isinstance(obj, np.ndarray):
                return obj.tolist()
            elif isinstance(obj, dict):
                return {k: convert(v) for k, v in obj.items()}
            elif isinstance(obj, list):
                return [convert(item) for item in obj]
            else:
                return obj

        converted_results = convert(results)

        flat_row = flatten_for_csv(converted_results)
        csv_out = output_dir / f"{Path(pdb_file).stem}_interface_results.csv"
        pd.DataFrame([flat_row]).to_csv(csv_out, index=False)
        
        print(f"\nFlat CSV saved to: {csv_out}")
        
        return 0
        
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
        return 1            

# For direct script execution
if __name__ == "__main__":
    
    exit(main())
