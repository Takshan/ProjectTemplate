```python
!pip install csfdock # Always install so that it gets updated.
```


```python
from csfdock.ar2a_v3 import *
```


```python
mol = ProjectStart()
```


```python
mol.SetFolders(".")  # sets curent working dir as Project Directory

```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">Project Base Directory: <span style="color: #800080; text-decoration-color: #800080">/home/lab09/SPACE/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff">Rahul-Iikwon</span>
</pre>




```python
mol.LoadReceptor("3eml_clean_receptor")

```


<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">Receptor: <span style="color: #800080; text-decoration-color: #800080; font-weight: bold">/home/lab09/SPACE/Rahul-Iikwon/DATA/</span><span style="color: #ff00ff; text-decoration-color: #ff00ff; font-weight: bold">3eml_clean_receptor.pdb</span>
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">
For<span style="color: #800000; text-decoration-color: #800000; font-weight: bold"> /home/lab09/SPACE/Rahul-Iikwon/DATA/3eml_clean_receptor.pdb</span>:
</pre>




<pre style="white-space:pre;overflow-x:auto;line-height:normal;font-family:Menlo,'DejaVu Sans Mono',consolas,'Courier New',monospace">┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┳━━━━━━━━━┓
┃<span style="color: #000080; text-decoration-color: #000080; font-weight: bold"> Record                                   </span>┃<span style="color: #000080; text-decoration-color: #000080; font-weight: bold">  Counts </span>┃
┡━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━╇━━━━━━━━━┩
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Chains:</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                                  </span>│       1 │
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Ligands:</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                                 </span>│       0 │
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Number of ligand atoms :</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                 </span>│       0 │
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Protein residues:</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                        </span>│     472 │
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Lipids molecules :</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                       </span>│       0 │
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Water molecules :</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                        </span>│       0 │
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Ions:</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                                    </span>│      96 │
│<span style="color: #7f7f7f; text-decoration-color: #7f7f7f"> </span><span style="color: #7fbf7f; text-decoration-color: #7fbf7f; font-weight: bold">Ion types :</span><span style="color: #7f7f7f; text-decoration-color: #7f7f7f">                              </span>│ {'HSD'} │
└──────────────────────────────────────────┴─────────┘
</pre>






    '/home/lab09/SPACE/Rahul-Iikwon/DATA/3eml_clean_receptor.pdb'




```python
mol.SimpleView()
```


