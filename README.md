- This is the README file for $(k,p)$-core maintenance code.

    The code can be run on Linux server with GCC 6.4.0 or later version.

    #### Setup and Run
    
    Enter the `kpcm` folder and run `make` command.

    If the command runs succesfully, a file `kpcore` will be generated.

    In the initial version of the program, the file `kpcore` contains the experiment : *Test the running time of the local algorithm under edge insertion*.

    Please use the following command to run the program, replace `<Dataset>` with name of the dataset you want to run tests on.
    
    ```shell
    ./kpcore <Dataset>
    ```
    
    The datasets we provided are listed below :
    
    ```
    Gnutella, Facebook, EmailEuAll, Gowalla, YouTube, Amazon, BerkStan
    ```
    
    Due to the 2GB storage constraints of Anonymous Github, we do not provide the datasets *LiveJournal* and *Orkut*.
    
    #### Experiments Configurations
    
    ##### Experiments Description
    
    You can run other experiments by replacing the statement `runtime_test_local_insertion(argv[1]);`  in the `main` function in the `main.cpp` file.
    
    You can find optional experiments in `experiment.h`. We list all functions and their corresponding experiments below.
    
    - `void correctness_test_local_insertion(std::string dataset)`
    
        *Test the correctness for the algorithm `Local` under edge insertions.*
    
    - `void correctness_test_local_deletion(std::string dataset)`
    
        *Test the correctness for the algorithm `Local` under edge deletions.*
    
    - `void runtime_test_local_insertion(std::string dataset)`
    
        *Test the runtime for the algorithm `Local` under edge insertions.*
    
    - `void runtime_test_local_deletion(std::string dataset)`
    
        *Test the runtime for the algorithm `Local` under edge deletions.*
    
    - `void runtime_test_local_woc_insertion(std::string dataset)`
    
        *Test the runtime for the algorithm `Local-woc` under edge insertions.*
    
    - `void runtime_test_local_woc_deletion(std::string dataset)`
    
        *Test the runtime for the algorithm `Local-woc` under edge deletions.*
    
    - `void vertices_statistics_local_insertion(std::string dataset)`
    
        *Test the vertices handled  for the algorithm `Local` under edge insertions, it will output the vertices maintained for both algorithm `pIncrease` and `pDecrease` with other statistic data , and also output the vertices with p-number changes in each insertion*.
    
    - `void vertices_statistics_local_deletion(std::string dataset)`
    
        *Test the vertices handled  for the algorithm `Local` under edge deletions, it will output the vertices maintained for both algorithm `pIncrease` and `pDecrease` with other statistic data , and also output the vertices with p-number changes in each deletion*.
    
    ##### An Example
    
    For example, you can replace the statement
    
    ```c++
    runtime_test_local_insertion(argv[1]);
    ```
    
    with 
    
    ```
    correctness_test_local_deletion(argv[1]);
    ```
    
    And then run `make` command to compile the file. Afterwards, using the command
    
    ```
    ./kpcore YouTube
    ```
    
    You can test the correctness of `Local` algorithm for edge deletion case on the dataset *Youtube*.
    
    #### Datasets & Generating Data
    
    ##### Datasets Provided
    
    The datasets we provided are listed below :
    
    ```
    Gnutella, Facebook, EmailEuAll, Gowalla, YouTube, Amazon, BerkStan
    ```
    
    Due to the 2GB storage constraints of Anonymous Github, we do not provide the datasets *LiveJournal* and *Orkut*.
    
    For the data used in scalability tests, we provide BerkStan's subgraphs of different scales.
    
    ```
    BerkStan_20E, BerkStan_40E, BerkStan_60E, BerkStan_80E
    ```
    
    The above four datasets are the subgraphs induced by 20% / 40% / 60% / 80% edges of the original BerkStan graph.
    
    ```
    BerkStan_20V, BerkStan_40V, BerkStan_60V, BerkStan_80V
    ```
    
    The above four datasets are the subgraphs induced by 20% / 40% / 60% / 80% vertices of the original BerkStan graph.
    
    ##### Generating Data
    
    We provide code for generating edges for insertion/deletion and code for generating subgraph induced by edges/vertices, the codes are in the directory `Tests`.
    
    - `generate_AddEdge.cpp`
    
        The code is used to generate edges for insertion. Run the following command to compile the code.
    
        ```shell
        g++ generate_AddEdge.cpp -o generate_AddEdge -std=c++11
        ```
    
        And the result program `generate_AddEdge` will be produced.
    
        Run the program using the following command
    
        ```
        ./generate_AddEdge <Dataset> <EdgeNum>
        ```
    
        Replace `<Dataset>` with dataset names provided, and replace `<EdgeNum>` with the number of edges you want to generate for insertion. The program will generate two files `_<Dataset>2.txt` and `_<Dataset>2_addGrph.txt` in the directory `datasets`, storing the graph without the edges for insertion, and the edges for insertion, respectively. The original files will be replaced by the new generated files.
    
    - `generate_DelEdge.cpp`
    
        The code is used to generate edges for deletion. Run the following command to compile the code.
    
        ```shell
        g++ generate_DelEdge.cpp -o generate_DelEdge -std=c++11
        ```
    
        And the result program `generate_DelEdge` will be produced.
    
        Run the program using the following command
    
        ```
        ./generate_DelEdge <Dataset> <EdgeNum>
        ```
    
        Replace `<Dataset>` with dataset names provided, and replace `<EdgeNum>` with the number of edges you want to generate for deletion. The program will generate two files `_<Dataset>full.txt` and `_<Dataset>full_delGrph.txt` in the directory `datasets`, storing the original graph, and the edges for deletion, respectively. The original files will be replaced by the new generated files.
    
    - `sample_graph_by_edge.cpp`
    
        The code is used to generate subgraphs induced by randomly seleted edges.  Run the following command to compile the code.
    
        ```shell
        g++ sample_graph_by_edge.cpp -o sample_graph_by_edge -std=c++11
        ```
    
        And the result program `sample_graph_by_edge` will be produced.
    
        Run the program using the following command
    
        ```
        ./sample_graph_by_edge <Dataset> <Proportion>
        ```
    
        Replace `<Dataset>` with dataset names provided, and replace `<EdgeNum>` with the a integer number in range $[0,100]$, which represents for the proportion of edges you want to choose (in persentage form). The program will generate a file`_<Dataset>_<Proportion>E.txt` , storing the graph induced by the randomly chosen edges.
    
    - `sample_graph_by_vertex.cpp`
    
        The code is used to generate subgraphs induced by randomly seleted vertices.  Run the following command to compile the code.
    
        ```shell
        g++ sample_graph_by_vertex.cpp -o sample_graph_by_vertex -std=c++11
        ```
    
        And the result program `sample_graph_by_vertex` will be produced.
    
        Run the program using the following command
    
        ```
        ./sample_graph_by_vertex <Dataset> <Proportion>
        ```
    
        Replace `<Dataset>` with dataset names provided, and replace `<EdgeNum>` with the a integer number in range $[0,100]$, which represents for the proportion of vertices you want to choose (in persentage form). The program will generate a file`_<Dataset>_<Proportion>V.txt` , storing the graph induced by the randomly chosen vertices.