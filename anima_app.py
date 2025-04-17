import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw
from st_social_media_links import SocialMediaIcons

from anima.models import AIkernel

with st.sidebar.container(border=False):
    st.sidebar.markdown(
        """
    # Anima
    - [Description](#desc)
    - [How to use](#htu)
    - [Anima predicitons](#pred)
    """,
        unsafe_allow_html=True,
    )

with st.sidebar.container(border=False):
    st.sidebar.markdown("""
    # About me
    You can contact me on one of the channels:
    """)
with st.sidebar.container(border=False):
    linkedin_image = "https://media.licdn.com/dms/image/v2/C4E03AQEYVPd4oVI9xg/profile-displayphoto-shrink_400_400/profile-displayphoto-shrink_400_400/0/1653659430741?e=1750291200&v=beta&t=_scrO_7DxXjfi4VfOBsZ1GWpYwE00b4wwZDrd9PkMvo"
    _, _, cent, _, _ = st.columns(5)
    with cent:
        st.image(linkedin_image, width=100)

with st.sidebar.container(border=False):
    social_media_links = [
        "https://www.linkedin.com/in/rpcarvs/",
        "mailto:rodrigo.carvalho.al@gmail.com",
    ]
    social_media_icons = SocialMediaIcons(
        social_media_links,
        colors=["#0A66C2", "#EF4026"],
    )
    social_media_icons.render()

with st.container(border=False):
    st.title("Anima")
    st.subheader("Description", divider=True, anchor="desc")
    st.markdown(
        """This page is a very simple visual and interactive representation
of the work done during my PhD. The idea was to develop
a robust AI-driven methodoly to boost the discovery of
novel organic-based electroactive materials for ion batteries.
"""
    )

    st.markdown("""
The process involved:

- Creating a large database of organic molecules comprising thousands of
entries, serving as a foundation for the framework.
- Designing a NLP-based machine learning model to predict the properties of interest,
completely by-passing time-intensive steps in materials research.
- Set a state-of-the-art framework to leverage the prediction models to efficiently
screen millions of candidate molecules, identifying high-performing battery materials.
""")
    st.markdown("""If you would like to know more about it, my work is completely described 
in my Doctoral Thesis, which can be access through the link below: """)

    st.link_button(
        "**My Thesis:** Organic Electrode Battery Materials: A Journey from Quantum Mechanics to Artificial Intelligence",
        "https://www.diva-portal.org/smash/record.jsf?pid=diva2%3A1687486&dswid=9698",
        use_container_width=True,
    )
    st.link_button(
        "**Paper:** Artificial intelligence driven in-silico discovery of novel organic lithium-ion battery cathodes",
        "http://dx.doi.org/10.1016/j.ensm.2021.10.029",
        use_container_width=True,
    )
    st.link_button(
        "**Paper:** An evolutionary-driven AI model discovering redox-stable organic electrode materials for alkali-ion batteries",
        "http://dx.doi.org/10.1016/j.ensm.2023.102865",
        use_container_width=True,
    )

    st.link_button(
        "Anima repository",
        "https://gitlab.com/rpcarvalho/anima",
        use_container_width=True,
    )
    st.caption("""This repo was organized on the final moments of my PhD,
               so do not expect high-quality or production-grade code 😄""")


with st.container(border=False):
    st.subheader("How to use", divider=True, anchor="htu")
    st.markdown("""Anima inputs must be molecules represented in the 
SMILES format. You can read more and maybe prepare some inputs using
the links below.""")

    col1, col2, col3 = st.columns(3)

    with col1:
        st.page_link(
            "https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System",
            label="**More about SMILES (Wikipedia)**",
        )

    with col2:
        st.page_link(
            "https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html",
            label="**SMILES generator and checker**",
        )
    with col3:
        st.page_link(
            "https://pubchem.ncbi.nlm.nih.gov//edit3/index.html",
            label="**SMILES sketcher from PubChem**",
        )

    st.markdown("""The SMILES string representing a given molecule
    (for example "CCCC=O") can be inserted on the designated field. After
    pressing 'Predict!', the input will be checked and adjusted to what the
    framework expects. If succeeds, the molecule will be ploted together with
    the predictions.""")

with st.container(border=False):
    st.subheader("Anima predicitons", divider=True, anchor="pred")

    st.markdown("""Add your SMILES and press 'Predict!'""")

    col1, col2 = st.columns([0.7, 0.3], vertical_alignment="bottom")
    with col1:
        molecule = st.text_input(
            "SMILES",
            value="BrC1=C2C(=O)N=Cc3c2n2C(O1)NC(=N)c2c3",
        )
    with col2:
        predict = st.button("Predict!")

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
            st.write("Ops! It seems you have intered an invalid SMILES")

        return result

    if predict and molecule:
        st.write(f"Molecule: {molecule}")

        # getting the kernel up and predictions
        results = run_model(molecule)
        new_molecule = results["smiles"][0]
        li_voltage = results["voltages"][0]
        oxidation = results["ox/red"][0][0]
        reduction = results["ox/red"][1][0]

        with st.container(border=True):
            st.write(f"**Molecule after preparation:** {new_molecule}")
            st.write(f"**Lithium insertion potential:** $${{\color{{green}}{round(float(li_voltage), 3)}}} \\textrm{{ V \\textit{{vs.}} Li/Li}}^+$$")
            st.write(f"**Molecule oxidation potential:** $${{\color{{red}}{round(float(oxidation), 3)}}} \\textrm{{ V (ref to vacuum)}}$$")
            st.write(f"**Molecule reduction potential:** $${{\color{{blue}}{round(float(reduction), 3)}}} \\textrm{{ V (ref to vacuum)}}$$")

            m = Chem.MolFromSmiles(new_molecule)
            fig = Draw.MolToImage(m)
            # plt.axis("off")
            fig
