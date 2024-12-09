Rscript ABCDmod.R 0 NSE 1 data3045.txt ABCDparam.txt

Los argumentos de la linea (y posibles valores):

    ABCDmod.R: nombre del archivo de código (el script principal del modelo)
    0: Modo de ejecución (si es ABCD por defecto -valor 0-, incluyendo nieve, etc -os lo explico al principio del script de código-) (posibles valores: 0 1 2 3)
    NSE: Función objetivo (posibles valores: MSE RMSE NSE R2 KGE)
    1: Si queremos guardar el modelo (puede ser 0 o 1)
    data3045.txt: Archivo de datos de entrada
    ABCDparam.txt: Archivo de parámetros de entrada. Son siete, aunque en función del tipo de ejecución coge unos u otros (ver explicación al princpio del archivo de código)

Como resultado se genera el archivo de ajuste con los siguientes datos (el úlitmo sería el valor a utilizar a minimizar -si es por ejemplo NSE aparece en negativo porque debería maximizarse-) (aunque podría haber dejado solamente un valor como me pedíais):

 

MSE,RMSE,PBIAS..,NSE,R2,KGE,f0
169.860563,13.033057,3.1,0.6098,0.742918,0.734632,-0.6098

 