# Computational High Energy Physics

Computational High Energy Physics is a field of study that focuses on using computational techniques and tools to understand and analyze the fundamental nature of the universe and the fundamental particles that make it up. This field is particularly relevant to the study of high energy physics, which involves the study of particles and phenomena that occur at extremely high energies, such as those found in particle accelerators and cosmic rays.

![](https://science.osti.gov/-/media/CR/Programs/High-Energy-Physics-Header-Image.jpg)


In Computational High Energy Physics, researchers use a range of techniques and tools, including numerical simulations, data analysis, and machine learning, to study and model the behavior of these fundamental particles and the forces that govern them. These techniques allow researchers to make predictions about the behavior of particles and forces under different conditions, and to test these predictions against experimental data.

Computational High Energy Physics plays a vital role in our understanding of the fundamental laws of the universe and the forces that shape it. It is also an important tool for the development of new technologies and advances in a wide range of fields, including medicine, materials science, and energy production.

---------


This repository contains projects and code relating to the field of Computational High Energy Physics, ranging from the simulation of **Lattice Quantum Chromodynamics** to Machine Learning and Quantum Machine Learning projects in the field of Experimental High Energy physics.






## Content

-------


### [Latttice Based Equation of State](https://github.com/MonitSharma/Computational-High-Energy-Physics/tree/main/Lattice%20Based%20Equation%20Of%20State) :

Lattice based equation of state (EOS) is a computational approach used to study the properties of matter under extreme conditions, such as those found in the cores of neutron stars or in the early universe. It is particularly relevant to the study of quantum chromodynamics (QCD), which is the theory that describes the interactions of quarks and gluons, the fundamental building blocks of matter.


![](https://itp.uni-frankfurt.de/~philipsen/grants/media/qcd_pd.png)


In lattice based EOS, researchers use numerical simulations to study the behavior of quarks and gluons on a lattice, a grid of points in space and time. This approach allows researchers to study the EOS of QCD under a wide range of conditions, including high temperatures and densities.

One of the main advantages of lattice based EOS is that it provides a first-principles approach to studying the EOS of QCD, meaning that it is based on the fundamental laws of the theory rather than on empirical models or assumptions. This makes it a powerful tool for understanding the properties of matter under extreme conditions and for testing the predictions of QCD against experimental data.

Lattice based EOS has been used to study a wide range of phenomena in QCD, including the phase structure of matter, the behavior of quark-gluon plasmas, and the properties of hadrons (particles made up of quarks and gluons). It is an active area of research with many open questions and challenges, and it continues to play a central role in our understanding of the fundamental nature of matter.


![](https://ars.els-cdn.com/content/image/1-s2.0-S0375947414003339-gr004.jpg)



-----------


This project simulates **Equation of State with Baryon Number, Electric charge** and **Strangeness Chemical Potential**. 

This is an implementation of the paper **"Lattice-based Equation of State at finite Baryon number, Electric charge and Strangeness chemical potential"** by J.Noronha-Hostler, P.Parotto, C. Ratti and J.M Stafford

Original paper : https://inspirehep.net/record/1720588

*This work should be cited whenever using the code.*

----

### [Equation of State with Critical Point QCD](https://github.com/MonitSharma/Computational-High-Energy-Physics/tree/main/EOS%20with%20Critical%20Point%20QCD):

The equation of state (EOS) is a relationship between the pressure, temperature, and density of a substance. In the context of quantum chromodynamics (QCD), the EOS describes the properties of quark-gluon plasma (QGP), a state of matter that is believed to have existed in the early universe and that can be produced in high-energy collisions in particle accelerators.

![](https://www.bnl.gov/riken/images/qcd-620px.jpg)

The EOS of QCD is of particular interest because it is expected to have a critical point, a singular point in the phase diagram at which the EOS exhibits a sudden change in behavior. This critical point is thought to be a sign of a phase transition, a change in the fundamental properties of matter, and it is believed to be associated with the deconfinement of quarks and gluons, the fundamental building blocks of matter.

There are many open questions and challenges in the study of the EOS of QCD with a critical point. Researchers are using a range of experimental and theoretical approaches to try to understand the properties of this phase transition and the nature of the QGP at high temperatures and densities. This includes the use of lattice based EOS, numerical simulations, and data analysis techniques, as well as the study of experimental data from particle accelerators.

Understanding the EOS of QCD with a critical point is of fundamental importance for our understanding of the fundamental nature of matter and the forces that govern it. It is also of practical importance for a wide range of applications, including the study of the properties of neutron stars and the behavior of matter under extreme conditions

![](https://www.usqcd.org/images/NewPhaseDiagram.png)

----- 

Lattice QCD simulations and finding critical point in QCD transitions

This is an updated version of the **BEST Collaboration** program producing an EoS matching lattice QCD at $\mu_B=0$, and containing a critical point in the $3D$ Ising universality class, in a parametrized form. It allows for different choices of constraints on the shape of the critical line, which reduce the number of parameters, as well as the inclusion of no critical point, corresponding to a Taylor expansion of lattice QCD result only.

Moreover, it allows for the choice of strangeness neutrality conditions: $n_S = 0$ && $n_Q = 0.4 n_B$ or the case of only baryon number: $\mu_S = \mu_Q = 0$

*When using this code, the following works should be cited: (https://inspirehep.net/literature/1672952) , (https://inspirehep.net/literature/1851774)*

--------

### [Addressing LHC problems via Machine Learning](https://github.com/MonitSharma/Computational-High-Energy-Physics/tree/main/Coursera%20Course%20Solutions):


![](https://cds.cern.ch/images/CERN-PHOTO-201802-030-10/file?size=large)


The Large Hadron Collider (LHC) is the largest data generation machine for the time being. It doesn’t produce the big data, the data is gigantic. Just one of the four experiments generates thousands gigabytes per second. The intensity of data flow is only going to be increased over the time. So the data processing techniques have to be quite sophisticated and unique.

In this online course students are introduced into the main concepts of the Physics behind those data flow so the main puzzles of the Universe Physicists are seeking answers for will be much more transparent. Of course the major stages of the data processing pipelines are scrutinized, and focus on the role of the Machine Learning techniques for such tasks as track pattern recognition, particle identification, online real-time processing (triggers) and search for very rare decays. The assignments of this course gives the opportunity to apply skills in the search for the New Physics using advanced data analysis techniques.

*The course can be found here: https://www.coursera.org/learn/hadron-collider-machine-learning*

--------


### [Deep Learning at CERN](https://github.com/MonitSharma/Computational-High-Energy-Physics/tree/main/Deep-Learning-in-High-Energy-Physics)

An Introduction to Deep Learning with Keras, focused on using it on data from CERN and classification and generation of jet images. data handling etc.


![](https://miro.medium.com/max/1070/1*bma4ZxCBH96atcM7ZQ0Wyg.jpeg)


At CERN, the European Organization for Nuclear Research, deep learning has been used in a number of research projects, including the analysis of data from particle accelerators and the search for new particles.

One example of the use of deep learning at CERN is in the analysis of data from the Large Hadron Collider (LHC), the world's largest and most powerful particle accelerator. The LHC generates a vast amount of data, and the use of deep learning techniques has allowed researchers to analyze this data more efficiently and effectively. For example, deep learning has been used to identify the signatures of new particles, such as the Higgs boson, in the data from the LHC.

![](https://knowledgetransfer.web.cern.ch/sites/default/files/styles/cern_image_gallery_large/public/images/competence/image/machine-learning-and-deep-learning.png?itok=HwAO5jH4)


Deep learning has also been used at CERN in the development of new particle detectors and in the optimization of accelerator beam performance. It has the potential to significantly improve the efficiency and precision of these systems and to enable new discoveries in high energy physics.

Overall, the use of deep learning at CERN has been an important tool for advancing our understanding of the fundamental nature of the universe and the fundamental particles that make it up. It is likely to continue to play a significant role in the research at CERN in the future.




-----


### [Quantum Computing Demo at CERN](https://github.com/MonitSharma/Computational-High-Energy-Physics/tree/main/Quantum%20Computing%20Demo%20%40%20CERN)


![](https://quanthep.eu/wp-content/uploads/2020/05/TheLOGO.svg)


Quantum computing is a new type of computing technology that makes use of the principles of quantum mechanics to perform calculations that are beyond the capabilities of classical computers. It has the potential to solve certain types of problems much more quickly and efficiently than classical computers, and it has attracted a great deal of interest from researchers and industries around the world.

![](https://scitechdaily.com/images/High-Energy-Physics-Quantum-Computing.jpg)


At CERN, the European Organization for Nuclear Research, quantum computing is being explored as a potential tool for a wide range of applications, including the simulation of particle collisions, the optimization of accelerator beam performance, and the analysis of large datasets.

One example of the use of quantum computing at CERN is in the simulation of particle collisions, which is a key part of the research at the organization. Quantum computers have the potential to simulate these collisions much more efficiently than classical computers, which could lead to a better understanding of the fundamental nature of the universe and the fundamental particles that make it up.

CERN is also exploring the use of quantum computing in the optimization of accelerator beam performance. Quantum computers have the potential to perform the calculations required to optimize the beam much more quickly than classical computers, which could lead to significant improvements in the efficiency and precision of these systems.

Overall, quantum computing has the potential to significantly advance the research at CERN and to enable new discoveries in high energy physics. It is an active area of research, and it is likely to continue to be an important focus of the organization in the future.

This repo contains materials from the Quantum Computing Demo given by Maria Schuld for CERN
This contains lecture slides and jupyter notebooks to get started with quantum computing.

*The even details can be found here: https://indico.cern.ch/event/893116/*



