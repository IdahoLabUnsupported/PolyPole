# PolyPole

*This project is not currently supported by Idaho National Laboratory as production grade software but is being made available to the public. Please use at your own risk.*

PolyPole is a numerical algorithm for the calculation of intra-granular fission gas release. In
particular, the algorithm solves the gas diffusion problem in a fuel grain in time-varying
conditions. The program has been extensively tested. PolyPole combines a high accuracy with a
high computational efficiency and is ideally suited for application in fuel performance codes.
The PolyPole algorithm is based on the analytic modal solution of the diffusion problem for
constant conditions, with the addition of polynomial corrective terms that embody the
information on the deviation from constant conditions. This semi-analytic concept is intended to
allow for significantly lower computational time than spatial discretization methods, such as
finite difference schemes. The PolyPole algorithm has been verified by comparing the results to
a finite difference reference solution for a large number of randomly generated operation
histories. Results demonstrated that (i) the accuracy of the PolyPole solution is high and
superior to other algorithms currently available, (ii) the computational time associated with
PolyPole is similar to other algorithms, and (iii) differently from other algorithms, the accuracy
of the PolyPole solution is highly consistent over the whole range of intra-granular fission gas
release.

Details of the PolyPole development and verification can be found in:

D. Pizzocri, C. Rabiti, L. Luzzi, T. Barani, P. Van Uffelen, G. Pastore. PolyPole-1: An accurate
numerical algorithm for intra-granular fission gas release, Journal of Nuclear Materials, 478,
333-342, 2016.

### Other Software
Idaho National Laboratory is a cutting edge research facility which is a constantly producing high quality research and software. Feel free to take a look at our other software and scientific offerings at:

[Primary Technology Offerings Page](https://www.inl.gov/inl-initiatives/technology-deployment)

[Supported Open Source Software](https://github.com/idaholab)

[Raw Experiment Open Source Software](https://github.com/IdahoLabResearch)

[Unsupported Open Source Software](https://github.com/IdahoLabCuttingBoard)

### License
Copyright 2017 Battelle Energy Alliance, LLC

Licensed under the Mozilla Public License version 2.0


