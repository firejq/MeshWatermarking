- [MeshWatermarking](#meshwatermarking)
  - [1. Dependence](#1-dependence)
  - [2. Watermarking Algorithm Implementation](#2-watermarking-algorithm-implementation)
    - [2.1. Embedding](#21-embedding)
    - [2.2. Extracting](#22-extracting)
    - [2.3. Correlation Checking](#23-correlation-checking)
  - [3. Watermarking Algorithm Improvement](#3-watermarking-algorithm-improvement)
    - [3.1. Watermarking sequence generation](#31-watermarking-sequence-generation)
    - [3.2. Embedding](#32-embedding)
    - [3.3. Extracting](#33-extracting)
  - [4. Experiment](#4-experiment)
    - [4.1. Mesh model](#41-mesh-model)
    - [4.2. Transparency](#42-transparency)
    - [4.3. Robustness](#43-robustness)
  - [5. Reference](#5-reference)
  - [6. License](#6-license)

# MeshWatermarking

This project implemented [a digital watermarking algorithm](http://graphicsinterface.org/proceedings/gi2001/gi2001-2/) for 3D Meshes based on spectral domain spectral coefficient perturbation. In addition, based on previous studies, the original algorithm is improved from watermark generation, watermark capacity and robustness against smooth attacks.

Experimental analysis shows that the improved algorithm can improve the watermark capacity and robustness against smooth attacks under the premise of maintaining good transparency.

## 1. Dependence

![image](http://img.cdn.firejq.com/jpg/2019/6/6/1eaed04b2f28cdb9250476a5cc2b5041.jpg)

- C++ 11
- Matlab R2012a
- Visual Studio 2015
- [OpenMesh 8.0 (static) 32-bit without apps](https://www.openmesh.org/media/Releases/8.0/OpenMesh-8.0-VS2015-32-Bit-no-apps.exe)
- [Eigen 3.1.0](https://bitbucket.org/eigen/eigen/get/3.1.0.zip)
- [Meshlab 2016](http://www.meshlab.net/#download)

## 2. Watermarking Algorithm Implementation

### 2.1. Embedding

The Watermarking is embed in the spectral coefficient:

1. Calculate the laplacian matrix of mesh model (n * n).
1. Execute the eigenvalue decomposition for the laplacian matrix and normalize the eigenvectors.
1. Project the coordinate of a vertex onto a normalized eigenvector so it will produce a mesh spectral coefficient of the vertex.
1. Embed the watermarking data into the spectral coefficient. Now consider watermarking i-th spectral coefficients of one of the spectral axes s:

    ![image](http://img.cdn.firejq.com/jpg/2019/6/6/0a4ed28514e93b4c8c9c509f4ab0f9b6.jpg)

    pi ∈ {−1, 1} is the pseudorandom number sequence (PRNS) generated from a known stego-key kw. And α(α > 0) is the modulation amplitude.

    Modulation processes for the other two spectral axes are identical.

1. Inverse-transforming the set of spectral coefficients back into the domain of vertex coordinates by using the following formula produces a watermarked polygonal mesh.

    ![image](http://img.cdn.firejq.com/jpg/2019/6/6/e51c70911d664dee3d8341601e16be38.jpg)

### 2.2. Extracting

![image](http://img.cdn.firejq.com/jpg/2019/6/6/12828a7188683770f53612063e5e8ef3.jpg)

### 2.3. Correlation Checking

![image](http://img.cdn.firejq.com/jpg/2019/5/14/9f50381f8cda95cfb18c80b65675eb81.jpg)

## 3. Watermarking Algorithm Improvement

### 3.1. Watermarking sequence generation

Generate the watermarking sequence by thresholding a image with realistic meaning (like logo):

![image](http://img.cdn.firejq.com/jpg/2019/6/6/40fa0d776b4681601bf36032313a4a3d.jpg)

### 3.2. Embedding

Approximately, smaller eigenvalues correspond to lower spatial frequencies, and larger eigenvalues correspond to higher spatial frequencies. Eigenvectors and spectral coefficients of the smaller eigenvalues represent global shape features, while eigenvectors and spectral coefficients of the larger eigenvalues represent local or detail shape features.

When the mesh model encounters a smooth transformation attack, the surface details of the model will be greatly damaged, and the overall shape can remain basically unchanged. Therefore, the spectral coefficient of the high-frequency component corresponding to the surface detail will also be greatly affected. However, the spectral coefficients of the low-frequency components corresponding to the overall appearance are less interfered.

In the original algorithm, when the watermark is embedded, the spectral coefficients of each vertex are directly embedded in the order of the vertex index, so the watermark information will be interfered to a large extent when subjected to a smooth attack.

If the watermark embedding can be preferentially embedded from the low frequency part of the grid, then when the number of model vertices is large enough, it can show better resistance to the smooth transformation.

![image](http://img.cdn.firejq.com/jpg/2019/6/6/4453f59f9f24cb02b153703f9d1062b4.jpg)

### 3.3. Extracting

Since the order has changed, this matrix is still not the original embedded watermark sequence. Therefore, similar to the logic at the time of embedding, it is necessary to construct a minimum heap based on the key value pairs of the vertex index and the eigenvalue according to the eigenvalue sequence of the Laplacian matrix of the original model, and take the current one from the minimum heap each time. The vertex index corresponding to the minimum eigenvalue is assigned, and the consecutive 3 bits corresponding to the index are assigned to the corresponding positions of the new watermark sequence matrix in the order of raster scanning. After the execution is completed, the watermark matrix is also converted into a one-dimensional watermark sequence in raster order.

![image](http://img.cdn.firejq.com/jpg/2019/6/6/f63a0fbd753665a046206a5bd26f3115.jpg)

## 4. Experiment

### 4.1. Mesh model

![image](http://img.cdn.firejq.com/jpg/2019/6/6/5f0e554d40797c604bfabbc2a9903d41.jpg)

### 4.2. Transparency

![image](http://img.cdn.firejq.com/jpg/2019/6/6/43d29166d1c6ee34b19914a0dbc93541.jpg)

### 4.3. Robustness

Taubin smoothing attacks:

![image](http://img.cdn.firejq.com/jpg/2019/6/6/4f095061d1b750f037a188b07fd2dfd5.jpg)

![image](http://img.cdn.firejq.com/jpg/2019/6/6/51b17684ca4c493d862767fba2d25c75.jpg)

## 5. Reference

[Ohbuchi R, Takahashi S, Miyazawa T, et al. Watermarking 3D polygonal meshes in the mesh spectral domain[C] Proceedings of Graphics interface.](http://graphicsinterface.org/proceedings/gi2001/gi2001-2/)

## 6. License

The MeshWatermarking is under the [GPL](https://github.com/firejq/MeshWatermarking/blob/master/LICENSE) License.