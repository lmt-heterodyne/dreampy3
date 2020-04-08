# DREAMPY3 Package
# Installation process

This installation is tested on Ubuntu 18.04 LTS.

1. Clone the repository
    ```
    git clone git@github.com:gopastro/dreampy3.git
    ```

2. Change into the dreampy3 directory
    ```
    cd dreampy3/
    ```

3. Create a virtual environment
    ```
    virtualenv -p python3 venv
    ```

4. Activate the virtual environment
    ```
    source venv/bin/activate
    ```

5. Install the required libraries.
    ```
    pip install -r requirements.txt
    ```

6. Setup ipython profile
   ```
   ipython profile create dreampy3
   ```

7. Launch dreampy3
   ```
   ipython --profile=dreampy3
   ```

