/** 
 * JanetStatus.java is an application that compiles and runs
 * under J2SE 5.0. It requires no other files.
 */


public class JanetStatus {
    public static void main(String[] args) {

    final int unrecoverableError = 8;
    final int errorInrecoveryProcess = 4;
    final int processing = 2;
    final int ready = 1;


    int status = 1;
    
 //  status = status | unrecoverableError;
    
    if (status == 1) {
         System.out.println("ready to receive requests");
    }
       
    }
}
