import streamlit as st
from st_social_media_links import SocialMediaIcons

from utils import show_output

with st.sidebar.container(border=False):
    st.sidebar.markdown(
        """
    # Anima
    - [Anima predicitons](#pred)
    - [How does it work?](#hwork)
    - [Further details](#desc)
    - [What are SMILES?](#sml)
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
    st.markdown(
        """<div style="text-align: justify;">
This page is a very simple visual and interactive representation
of the work done during my PhD. Add your SMILES molecule and press 
Predict!
</div>""",
        unsafe_allow_html=True,
    )
    st.subheader("Predicitons", divider=True, anchor="pred")
    # O=CC(=O)C(=O)N=C1C(C(=O)C=O)=CN=C1Br
    # N=c1sc2c(C(=O)C(=O)C=O)nc(=O)c=2c1=O
    st.markdown("""Add your SMILES and press **'Predict!**'""")
    st.markdown("""Some suggestions:

- O=CC(=O)C(=O)N=C1C(C(=O)C=O)=CN=C1Br
- N=c1sc2c(C(=O)C(=O)C=O)nc(=O)c=2c1=O
- BrC1=C2C(=O)N=Cc3c2n2C(O1)NC(=N)c2c3

 The max SMILES length is 54 elements.
""")

    col1, col2 = st.columns([0.7, 0.3], vertical_alignment="bottom")
    with col1:
        molecule = st.text_input(
            "**SMILES**",
            value="C1=CC(=CC=C1C(=O)O)C(=O)O",
        )
    with col2:
        predict = st.button("**Predict!**", type="primary", use_container_width=True)

    if predict and molecule:
        show_output(molecule)
    else:
        show_output("C1=CC(=CC=C1C(=O)O)C(=O)O")


with st.container(border=False):
    st.subheader("How does it work?", divider=True, anchor="hwork")
    st.markdown(
        """<div style="text-align: justify;">
The NLP Neural Network first breaks the SMILES into smaller pieces (tokens). These 
pieces, representing individual atom species, special bonds and others, are uniquely 
identified. Thus, a carbon is uniquely identified in this scheme, although a normal 
carbon is different than an aromatic carbon. These are then transformed into vector 
embeddings, further learned during training. This gives the model an extra level of 
abstraction to understand "what is a carbon?", "what is an oxygen?" or even "what is 
a double bond?" in the context of redox potentials. The embeddings are then processed 
through a parallel stack of recurrent (GRU) and dense layers. These are then concatenated 
and further decoded in a dense layer. This process is summarized in the figure below.
</div>""",
        unsafe_allow_html=True,
    )
    st.write("")
    st.image(
        image="./static/panel_smiles_nn.png",
        caption="Schematics of how the NLP-based SMILES Neural Net works. Figure extracted from [this paper](https://doi.org/10.1016/j.ensm.2021.10.029)",
        use_container_width=True,
    )
    st.write("")

    st.markdown(
        """<div style="text-align: justify;">
The next figure shows how the final kernel model works. Combined with a simpler regression 
model (linear model), the NLP neural net (Neural Model) is capable of predicing the redox 
potential (oxidation and reduction) and the lithiation (or the Lithium insertion) potential for 
the target molecule. This kernel can, for example, be easily act as an AI-Agent tool providing 
novel functionalities in the materials discovery process. The AI-Agent can autonomously search 
for new SMILES, potentially discovering new materials.
</div>""",
        unsafe_allow_html=True,
    )
    st.write("")
    st.image(
        image="./static/new_scheme.jpg",
        caption="Illustrating the final workflow. Figure extracted from [this paper](https://doi.org/10.1016/j.ensm.2021.10.029)",
        use_container_width=True,
    )
    st.write("")

with st.container(border=False):
    st.subheader("Futher details", divider=True, anchor="desc")
    st.markdown(
        """<div style="text-align: justify;">
The idea behind this project was to develop
a robust AI-driven methodoly to boost the discovery of
novel organic-based electroactive materials for ion batteries.
</div>""",
        unsafe_allow_html=True,
    )

    st.markdown(
        """<div style="text-align: justify;">
The process involved:

- Creating a large database of organic molecules comprising thousands of
entries, serving as a foundation for the framework.
- Designing a NLP-based machine learning model to predict the properties of interest,
completely by-passing time-intensive steps in materials research.
- Set a state-of-the-art framework to leverage the prediction models to efficiently
screen millions of candidate molecules, identifying high-performing battery materials.
</div>""",
        unsafe_allow_html=True,
    )
    st.markdown(
        """<div style="text-align: justify;">If you would like to know more about it, my work is completely described 
in my Doctoral Thesis, which can be access through the link below: </div>""",
        unsafe_allow_html=True,
    )

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
        "https://github.com/rpcarvs/anima",
        use_container_width=True,
        type="primary",
    )
    st.caption("""This repo was organized on the final moments of my PhD,
               so do not expect production-grade code 😄. In fact, it is quite messy.""")


with st.container(border=False):
    st.subheader("What are SMILES?", divider=True, anchor="sml")
    st.markdown(
        """<div style="text-align: justify;">Anima inputs must be molecules represented in the 
SMILES format. You can read more and maybe prepare some inputs using
the links below.</div>""",
        unsafe_allow_html=True,
    )

    col1, col2, col3 = st.columns(3)

    with col1:
        st.link_button(
            url="https://en.wikipedia.org/wiki/Simplified_Molecular_Input_Line_Entry_System",
            label="**More about SMILES (Wikipedia)**",
        )

    with col2:
        st.link_button(
            url="https://www.cheminfo.org/flavor/malaria/Utilities/SMILES_generator___checker/index.html",
            label="**SMILES generator and checker**",
        )
    with col3:
        st.link_button(
            url="https://pubchem.ncbi.nlm.nih.gov//edit3/index.html",
            label="**SMILES sketcher from PubChem**",
        )

    st.markdown(
        """<div style="text-align: justify;">The SMILES string representing a given molecule
    (for example "CCCC=O") can be inserted on the designated field. After
    pressing 'Predict!', the input will be checked and adjusted to what the
    framework expects. If successful, the molecule will be ploted together with
    the predictions.</div>""",
        unsafe_allow_html=True,
    )
