import streamlit as st
from anima.models import AIkernel
from rdkit import Chem
from rdkit.Chem import Draw


@st.cache_data
def run_model(smiles):
    try:
        result = AIkernel(
            [smiles],
            return_redox=True,
            return_smiles=True,
        )
    except Exception:
        result = {}
        st.write("Ops! It seems you have an invalid SMILES")

    return result


def show_output(molecule):
    st.write(f"Molecule: {molecule}")

    # getting the kernel up and predictions
    results = run_model(molecule)
    new_molecule = results["smiles"][0]  # type: ignore
    li_voltage = results["voltages"][0]  # type: ignore
    oxidation = results["ox/red"][0][0]  # type: ignore
    reduction = results["ox/red"][1][0]  # type: ignore

    with st.container(border=True):
        st.write(f"**Molecule after preparation:** {new_molecule}")
        st.write(
            f"**Lithium insertion potential:** $${{\color{{green}}{round(float(li_voltage), 3)}}} \\textrm{{ V \\textit{{vs.}} Li/Li}}^+$$"
        )
        st.write(
            f"**Molecule oxidation potential:** $${{\color{{red}}{round(float(oxidation), 3)}}} \\textrm{{ V (ref to vacuum)}}$$"
        )
        st.write(
            f"**Molecule reduction potential:** $${{\color{{blue}}{round(float(reduction), 3)}}} \\textrm{{ V (ref to vacuum)}}$$"
        )

        m = Chem.MolFromSmiles(new_molecule)
        fig = Draw.MolToImage(m)
        # plt.axis("off")
        st.image(fig)
