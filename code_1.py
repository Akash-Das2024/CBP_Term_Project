from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.SASA import ShrakeRupley
import matplotlib.pyplot as plt
import os


class ProteinSASA:
  """
  Class to calculate and analyze Solvent Accessible Surface Area (SASA) of a protein at different levels.
  """

  def __init__(self, pdb_file, probe_radius=1.4):
    """
    Initializes the ProteinSASA object with PDB file path and probe radius.

    Args:
        pdb_file (str): Path to the PDB file containing the protein structure.
        probe_radius (float, optional): Radius of the probe sphere (in Angstrom). Defaults to 1.4.
    """

    self.pdb_file = pdb_file
    self.probe_radius = probe_radius
    self.sasa_dict = {}

  def calculate_sasa(self):
    """
    Calculates the SASA for the protein at all levels (atoms, residues, chains, total).
    """

    filename = os.path.splitext(os.path.basename(self.pdb_file))[0]
    p = PDBParser(QUIET=1)
    struct = p.get_structure(filename, self.pdb_file)
    sr = ShrakeRupley()
    for level in ("A", "R", "C", "M", "S"):
      sr.compute(struct, level)
      self.sasa_dict[level] = round(struct.sasa, 2)

  def get_total_sasa(self, level):
    """
    Returns the total SASA for the protein at a specified level.

    Args:
        level (str): Level for which total SASA is requested (e.g., "A" for atoms, "R" for residues, "C" for chains, "total").

    Returns:
        float: Total SASA for the protein at the specified level.
    """

    if not self.sasa_dict:
      self.calculate_sasa()

    if level not in self.sasa_dict:
      raise ValueError(f"Invalid level '{level}'. Supported levels: A, R, C, M, S")

    return self.sasa_dict[level]

  def plot_sasa(self, level, save_plot=True):
    """
    Plots the SASA for the protein at a specified level.

    Args:
        level (str): Level for which SASA is plotted (e.g., "A" for atoms, "R" for residues, "C" for chains).
    """

    if not self.sasa_dict:
      self.calculate_sasa()

    if level not in self.sasa_dict:
      raise ValueError(f"Invalid level '{level}'. Supported levels: A, R, C")

    if level == "A":
      data_x, data_y = zip(*self.sasa_dict[level].items())
      plt.bar(data_x, data_y)
      plt.xlabel("Atom ID")
    elif level == "R":
      data_x, data_y = zip(*self.sasa_dict[level].items())
      plt.bar(data_x, data_y)
      plt.xlabel("Residue ID")
    elif level == "C":
      data_x, data_y = zip(*self.sasa_dict[level].items())
      plt.bar(data_x, data_y)
      plt.xlabel("Chain ID")
    plt.ylabel("SASA (Angstrom^2)")
    plt.title(f"{level} Level SASA")
    if save_plot:
      # Create Results folder if it doesn't exist
      results_folder = "Results"
      os.makedirs(results_folder, exist_ok=True)  # Create folder if it doesn't exist

      # Extract filename from PDB file path
      filename = os.path.splitext(os.path.basename(self.pdb_file))[0]
      plot_filepath = os.path.join(results_folder, f"{filename}_{level}_sasa.png")

      plt.savefig(plot_filepath)  # Save plot with filename based on level
      print(f"Plot saved to: {plot_filepath}")
    else:
      plt.show()
      plt.close()  # Close the plot window

  def print_results(self):
    """
    Prints the total SASA for the protein at all levels (atoms, residues, chains, total).
    """

    if not self.sasa_dict:
      self.calculate_sasa()

    for level in ("A", "R", "C", "M", "S"):
      total_sasa = self.get_total_sasa(level)
      level_name = {
        "A": "Atoms",
        "R": "Residues",
        "C": "Chains",
        "M": "Model",
        "S": "Structure"
      }[level]
      print(f"{level_name} SASA: {total_sasa:.2f} Å²")
    print("-"*20)


def main(pdb_folder):
  """
  Analyzes SASA for all PDB files in a specified folder.

  Args:
      pdb_folder (str): Path to the folder containing PDB files.
  """

  for filename in os.listdir(pdb_folder):
    if filename.endswith(".pdb"):
      pdb_filepath = os.path.join(pdb_folder, filename)
      protein = ProteinSASA(pdb_filepath)
      protein.calculate_sasa()
      print(f"Protein: {filename}")
      protein.print_results()
      for level in ("A", "R", "C", "M", "S"):
        protein.plot_residue_sasa(level)
      print("-"*20)

if __name__ == "__main__":
  pdb_folder = "PDB"
  main(pdb_folder=pdb_folder)
